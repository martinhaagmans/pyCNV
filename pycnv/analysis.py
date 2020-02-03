#!/usr/bin/env python
import os
import sys

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from ngsscriptlibrary import SexForCNV

from .plots import SeriePlots
from .plots import SamplePlots
from .databases import Databases

SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))


def get_annot_bedlocation(capture, pipelinedir=None):
    """Find annotated BEDfile in pipeline home. Return a filelocation."""
    bed_annot = os.path.join(pipelinedir, 'captures',
                             '{}_target.annotated'.format(capture))
    if not os.path.isfile(bed_annot):
        raise IOError('{} does not exist.'.format(bed_annot))
    return bed_annot


def annotbed_to_df(bedfile):
    """Read BEDfile with annotation info and return a dataframe"""
    dfannot = pd.read_csv(bedfile, header=None, sep='\t')
    try:
        dfannot.drop('strand', inplace=True, axis=1)
    except KeyError:
        pass
    if len(dfannot.columns) == 4:
        dfannot.columns = ['chromosome', 'start', 'end', 'Gen']
    elif len(dfannot.columns) == 5:
        dfannot.columns = ['chromosome', 'start', 'end', 'Gen', 'strand']

    dfannot['start'] = dfannot['start'] + 1
    dfannot.set_index(dfannot.apply(
        lambda x: '{}_{}_{}'.format(x['chromosome'],
                                    x['start'],
                                    x['end']), axis=1),
                      inplace=True)
    dfannot.drop(['chromosome', 'start', 'end'],
                 inplace=True, axis=1)
    dfannot.sort_index(inplace=True)
    return dfannot


def create_capture_database(capture, configfile=None):
    """Create 1 database with 2 tables:
     - DoC table for coverage data
     - Annotation table with gene-target info.
     """
    if configfile is None:
        configfile = os.path.join(SCRIPTDIR, 'config.py')
    config = get_config_dict(configfile)
    annotbed = get_annot_bedlocation(capture, pipelinedir=config['pipelinedir'])
    dfannot = annotbed_to_df(annotbed)
    DB = Databases(capture)
    DB.create_annot_table(dfannot)
    DB.create_doc_table()
    return


def underscore_targets(targets):
    targets = targets.str.replace(':', '_')
    targets = targets.str.replace('-', '_')
    return targets


def remove_underscore_targets(targets):
    targets = targets.str.replace('_', ':', 1)
    targets = targets.str.replace('_', '-', 1)
    return targets


def normalize_df(df, correctsex=True):
    """Get ratio target_coverage/average_coverage for each target for each
    patient. Return dataframe.
    """
    if type(df) == pd.DataFrame:
        patmean = pd.DataFrame(df.mean(axis=0)).transpose()
        dfnorm = df.div(patmean.iloc[0], axis='columns')
    elif type(df) == pd.Series:
        dfnorm = df.div(df.mean())

    return dfnorm.transpose()


def get_target_info(df):
    """Get mean and sd per target. Return dataframe."""
    dfnorm = normalize_df(df)
    targetmean = pd.DataFrame(dfnorm.mean(axis=0))
    targetmean.columns = ['Mean']
    targetstd = pd.DataFrame(dfnorm.std(axis=0))
    targetstd.columns = ['Std']
    return (targetmean, targetstd)


def get_zscore_df(df):
    """Normalizes dataframe and determines Z-scores."""
    dfnorm = normalize_df(df)
    targetmean, targetstd = get_target_info(df)
    targetstd.replace(0, 1, inplace=True)
    zscore = dfnorm.subtract(targetmean.transpose().iloc[0], axis=1)
    zscore = zscore.div(targetstd.transpose().iloc[0], axis=1).transpose()
    return zscore.round(1)


def correct_males(df, patientinfo_db, newdir, cutoff=0.85):
    """Calculate ratio X-chromosome coverage and Autosomal coverage and
    double the coverage on the X-chromosome if ratio < cutoff.
    Return a corrected dataframe
    """
    sex_report = '{}/geschat_geslacht.txt'.format(newdir)

    df = df.transpose()
    if not df.index.str.contains('^chrX', regex=True, na=False).any():
        return df.transpose()

    sex_unknown = list()
    males = list()
    
    S = SexForCNV(patientinfo_db)

    sex = df[df.index.str.contains('^chrX', regex=True, na=False)]
    auto = df[~df.index.str.contains('^chrX', regex=True, na=False)]

    for col in df.columns:
        is_male = None
        serie, sample = col

        if not S.sample_in_db(serie, sample):
            sex_unknown.append(col)
            continue
        try:
            is_male = S.sample_is_male(serie, sample)
        except ValueError:
            sex_unknown.append(col)
            continue

        if is_male:
            males.append(col)

    with open(sex_report, 'w') as f:      
        for col in sex_unknown:
            serie, sample = col
            f.write('{}\t'.format(sample))
            if (sex[sex[col] > 100][col].mean() / auto[auto[col] > 100][col].mean()) < cutoff:
                males.append(col)
                f.write('M\n')
            else:
                f.write('V\n')

    for col in males:
        new = 2 * df[df.index.str.contains('^chrX', regex=True,
                                           na=False)][col].values

        df.at[df.index.str.contains('^chrX', regex=True, na=False), col] = new
    return df.transpose()


def get_bad_regions(df, meancutoff=0.2, stdcutoff=0.15):
    """Define less/non callable regions by mean and std. Return dataframe."""
    excluded_std = df[df['Std'] >= stdcutoff]['Std']
    excluded_mean = df[df['Mean'] <= meancutoff]['Mean']
    excluded_regions = pd.concat([excluded_mean, excluded_std],
                                 axis=1, sort=True).fillna('OK')
    return(excluded_regions)


def percentage_callable(df, badregions):
    df.columns = ['target', 'gen']
    df['chr'], df['start'], df['end'] = zip(*df['target'].apply(lambda x: x.split('_')))
    df[['start', 'end']] = df[['start', 'end']].apply(pd.to_numeric)
    df['size'] = df['end'] - df['start']

    badregions['chr'], badregions['start'], badregions['end'] = zip(*badregions['index'].apply(lambda x: x.split('_')))
    badregions[['start', 'end']] = badregions[['start', 'end']].apply(pd.to_numeric)
    badregions['size'] = badregions['end'] - badregions['start']
    return (badregions['size'].sum() / df['size'].sum())


def create_database(args, configfile=None):
    if args.poscontroles or args.ingestuurd:
        DB = Databases('', configfile=configfile)
        if args.poscontroles:
            DB.add_poscontrols(args.poscontroles)
            sys.exit()
        if args.ingestuurd:
            DB.add_badsamples(args.ingestuurd)
            sys.exit()
    if args.capture:
        create_capture_database(args.capture, configfile=configfile)


def add_docfile(docfile, capture, serie, sample):
    D = Databases(capture)
    data = D.parse_docfile(docfile)
    D.add_data_to_db(sample, serie, data)


def drop_badsamples(dfdrop, badsampleIDs):
    try:
        drop = [i for i in badsampleIDs if i in dfdrop.index.droplevel(0)]
    except (AttributeError, ValueError):
        [dfdrop.drop(i, inplace=True) for i in badsampleIDs if i in dfdrop.index]
    else:
        dfdrop = dfdrop.drop(drop, level=1)
    return dfdrop


def drop_poscons(df, posconIDs):
    try:
        drop = [i for i in posconIDs if i in df.index.droplevel(0)]
    except (AttributeError, ValueError):
        [df.drop(i, inplace=True) for i in posconIDs if i in df.index]
        dfclean = df
    else:
        dfclean = df.drop(drop, level=1)
    return dfclean


def get_samples_for_serie(df, serie):
    samples = list(df.loc[serie].index)
    return samples


def get_poscondata(df, posconsamples):
    dfall = df.copy()
    try:
        df_poscon = dfall.loc[posconsamples]
    except KeyError:
        return pd.DataFrame()

    return df_poscon


def seperate_data(df, new, sample=None):
    """Seperate dataframe in archive and to be analyzed (new) data.
    Return 2 dataframes.
    """
    if sample is not None:
        df_new = df.loc[(new, sample)]
        df.drop([sample], level=1, inplace=True)
    else:
        df_new = df.loc[new]
        df.drop([new], level=0, inplace=True)
    return (df_new, df)


def get_intervals_for_gene(df_annot, gene):
    """Get intervals for gene from annotdf. Returns (list/index)"""
    return df_annot[df_annot['gen'].str.contains('^{}$'.format(gene))].index


def filter_data_by_intervals(data, intervals):
    """Filter dataframe by intervals. Return dataframe."""
    return data[data.index.isin(intervals)]


def sort_by_interval(df):
    df['target'] = df.index
    df['chr'], df['start'], df['end'] = zip(*df['target'].apply(lambda x: x.split('_')))
    df[['start', 'end']] = df[['start', 'end']].apply(pd.to_numeric)
    df.sort_values('start', inplace=True)
    df.drop(['target', 'chr', 'start', 'end'], axis=1, inplace=True)
    return df


def create_dirs(output, capture, serie, outdir):
    if output is not None:
        newdir = output
    elif output is None:
        newdir = os.path.join(outdir, capture, serie)
    newperpatcalls = os.path.join(newdir, 'Calls')
    newperpatqc = os.path.join(newdir, 'QC')
    PATHS = [newdir, newperpatcalls, newperpatqc]
    [os.makedirs(path) for path in PATHS if not os.path.exists(path)]
    return newdir


def get_config_dict(configfile):
    config = {}
    exec(open(configfile).read(), config)
    return config


def get_gene_list_from_file(genefile):
    with open(genefile, 'r') as f:
        return [line.strip() for line in f]


def serie_qc(df, capture, serie, outdir, poscons, badsamples):
    df = drop_badsamples(df, badsamples)
    df_new, df_archive = seperate_data(df, serie)
    df_archive = drop_poscons(df_archive, poscons)
    archive_nr_of_samples = len(df_archive.index)
    archive_targetmean, archive_targetstd = get_target_info(df_archive.transpose())
    archive_info = archive_targetmean.join(archive_targetstd)

    serie_targetmean, serie_targetstd = get_target_info(df_new.transpose())
    serie_info = serie_targetmean.join(serie_targetstd)
    pdf = PdfPages('{}/{}.pdf'.format(outdir, serie))
    Plotter = SeriePlots(capture, serie, pdf)
    Plotter.plot_qc_serie(serie_info,
                          archive_info,
                          archive_nr_of_samples
                          )
    pdf.close()
    print('{} QC done'.format(serie))


def sample_qc():
    pass


def write_archive_file(newdir, samples):
    with open('{}/archive.txt'.format(newdir), 'w') as f:
        for sample in samples:
            f.write(f'{sample}\n')


def write_excluded_file(newdir, badregions, empiricalfragments, perc_callable):
    if not badregions.empty:
        badregions.index = remove_underscore_targets(badregions.index)
        badregions.to_csv('{}/excluded.txt'.format(newdir), sep='\t')
        badregions.index = underscore_targets(badregions.index)

    with open('{}/excluded.txt'.format(newdir), 'a') as f:
        for frags in empiricalfragments:
            for frag in frags.split(', '):
                frag.replace('_', ':', 0)
                frag.replace('_', '-', 0)
                f.write(f'{frag}\tEmpirisch bepaald\n')

        if not badregions.empty:
            f.write('\n{0:0.1f} % is niet callable'.format(perc_callable))
        if badregions.empty:
            f.write('Geen.')


def analyse(capture, serie, docfile=None, sample=None, outdir=None,
            reportgenes=None, addonly=False, configfile=None, delete=False):
    if docfile:
        add_docfile(docfile, capture, serie, sample)
        if addonly:
            sys.exit()

    config = get_config_dict(os.path.join(SCRIPTDIR, 'config.py'))
    
    if delete:
        Databases(capture).delete_serie(serie)
        sys.exit()
    
    if not outdir:
        outdir = config['outputdir']
        newdir = create_dirs(None, capture, serie, outdir)
    elif outdir:
        newdir = create_dirs(outdir, capture, serie, outdir)    
    
    QD = Databases(capture)
    df = QD.get_archive()
    df = correct_males(df, config['patientinfo'], newdir)

    poscon_dict = QD.get_positive_controls_dict()
    poscon_ids = list(poscon_dict.keys())
    badsamples = QD.get_bad_samples()
    df_annot = QD.get_annot()
    empiricalfragments = QD.get_regions_to_exclude()



    serie_qc(df, capture, serie, newdir, poscon_ids, badsamples)

    serie_samples = get_samples_for_serie(df, serie)
    badsamples_archive = [_ for _ in badsamples if _ not in serie_samples]

    df = drop_badsamples(df, badsamples_archive)
    df.index = df.index.droplevel(0)
    df_poscons = df.reindex(poscon_ids)

    try:
        df.drop(poscon_ids, inplace=True)
    except KeyError as e:
        df_poscons = pd.DataFrame()

    df_clean = df.copy()
    df_clean = drop_badsamples(df_clean, badsamples)
    df_clean = df_clean.transpose()
    df = df.transpose()
    df_poscons = df_poscons.transpose()

    target_mean, target_std = get_target_info(df_clean)
    target_info = target_mean.join(target_std)

    all_samples = list(df_clean.columns)
    write_archive_file(newdir, all_samples)

    badregions = get_bad_regions(target_info)
    badregions = badregions.join(df_annot)
    if badregions.empty:
        perc_callable = 100
    elif not badregions.empty:
        perc_callable = percentage_callable(df_annot.reset_index(),
                                            badregions.reset_index())
        perc_callable = perc_callable * 100
    write_excluded_file(newdir, badregions, empiricalfragments, perc_callable)

    if not df_poscons.empty:
        df_poscon_zscore_list = list()
        for posconid in poscon_ids:
            df_temp = df.join(df_poscons[posconid])
            zscores_temp = get_zscore_df(df_temp)
            df_poscon_zscore_list.append(zscores_temp[posconid])

        df_zscores_poscons = pd.concat(df_poscon_zscore_list, axis=1)
        df_zscores_samples = get_zscore_df(df)
        df_zscores = df_zscores_samples.join(df_zscores_poscons)

    elif df_poscons.empty:
        df_zscores = get_zscore_df(df)

    for i, sample in enumerate(serie_samples):
        samples_to_drop = [_ for _ in badsamples if _ != sample]
        df_temp = drop_badsamples(df, samples_to_drop)
        df_norm = normalize_df(df_temp).transpose()
        sampleplotdf = target_mean.join(df_norm.join(df_annot))
        SaPlotter = SamplePlots(sample, newdir, badsamples=badsamples)
        SaPlotter.sample_qc(sampleplotdf, serie)
        calls = df_zscores[(df_zscores[sample] < -3) | (df_zscores[sample] > 3)]

        if len(calls.index) != 0:
            calls = calls.transpose()
            calls_per_target = {target: '{}/{}'.format(len(calls[(calls[target] > 3) | (calls[target] < -3)][target]),
                                                       len(calls[target]))
                                for target in calls}

            calls = calls.transpose()
            calls = df_annot.join(calls, how='right')
            genes = calls['gen'].dropna().unique()
            calls = calls[[sample, 'gen']].join(pd.DataFrame.from_dict(
                calls_per_target, orient='index')
                ).join(badregions[['Mean', 'Std']]).fillna('OK')
            calls.columns = ['Z-score', 'Gen', 'Freq', 'Mean', 'Std']
            calls.index = remove_underscore_targets(calls.index)
            calls.to_csv('{}/Calls/{}.txt'.format(newdir, sample), sep='\t')
            samplepdf = PdfPages('{}/Calls/{}.pdf'.format(newdir, sample))

            for gene in genes:
                if reportgenes is not None and gene not in reportgenes:
                    continue
                elif gene == 'nan':
                    continue
                try:
                    intervalstoplot = get_intervals_for_gene(df_annot, gene)
                except ValueError as e:
                    print(e, gene)
                    continue
                datatoplot = filter_data_by_intervals(df_zscores, intervalstoplot)
                try:
                    datatoplot = sort_by_interval(datatoplot.copy())
                except ValueError as e:
                    print(e, gene)
                    continue

                targetinfo_toplot = filter_data_by_intervals(target_info,
                                                             intervalstoplot)
                targetinfo_toplot = sort_by_interval(targetinfo_toplot.copy())

                SaPlotter.plot_cnv_calls(datatoplot, gene, samplepdf,
                                         targetinfo_toplot, serie_samples, poscon_dict)

            samplepdf.close()
        print('{} ({} of {}) done'.format(sample, i + 1, len(serie_samples)))
