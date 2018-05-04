#!/usr/bin/env python
import os
import sys

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .databases import Databases
from .plots import SeriePlots
from .plots import SamplePlots


def doc_to_df(docfile):
    """Read GATK's DoC file and return a df."""
    df = pd.read_csv(docfile, index_col=0, sep='\t')
    df = df.filter(regex='_mean_cvg')
    df.columns = [col.split("_")[0] for col in df]
    df.index = underscore_targets(df.index)
    df = df.transpose()
    df.sort_index(inplace=True)
    return df


def annotbed_to_df(bedfile):
    """Read BEDfile with annotation info and return a dataframe"""
    dfannot = pd.read_csv(bedfile, header=None, sep='\t')
    try:
        dfannot.drop('strand', inplace=True, axis=1)
    except ValueError:
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


def correct_males(df, cutoff=0.85):
    """Calculate ratio X-chromosome coverage and Autosomal coverage and
    double the coverage on the X-chromosome if ratio < cutoff.
    Return a corrected dataframe
    """
    df = df.transpose()
    sex = df[df.index.str.contains('^chrX', regex=True, na=False)]
    auto = df[~df.index.str.contains('^chrX', regex=True, na=False)]
    males = [x for x in df.columns
             if (sex[sex[x] > 100][x].mean() /
                 auto[auto[x] > 100][x].mean()) < cutoff]

    for p in males:
        new = 2 * df[df.index.str.contains('^chrX', regex=True,
                                           na=False)][p].values

        df.at[df.index.str.contains('^chrX', regex=True, na=False), p] = new
    return df.transpose()


def get_bad_regions(df, meancutoff=0.2, stdcutoff=0.15):
    """Define less/non callable regions by mean and std. Return dataframe."""
    targetmean, targetstd = get_target_info(df)
    excluded_std = targetstd[targetstd['Std'] >= stdcutoff]
    excluded_mean = targetmean[targetmean['Mean'] <= meancutoff]
    excluded_regions = pd.concat([excluded_mean, excluded_std],
                                 axis=1).fillna('OK')
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


def create_database(args):
    if args.poscontroles or args.ingestuurd:
        DB = Databases('')
        if args.poscontroles:
            DB.add_poscontrols(args.poscontroles)
        if args.ingestuurd:
            DB.add_badsamples(args.ingestuurd)
    if args.capture:
        create_capture_database(args.capture)


def create_capture_database(capture):
    """Create 1 database with 2 tables:
     - DoC table for coverage data
     - Annotation table with gene-target info.
     """
    configfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'config.py')
    config = get_config_dict(configfile)
    annotbed = get_annot_bedlocation(capture,
                                     pipelinedir=config['pipelinedir'])
    dfannot = annotbed_to_df(annotbed)
    DB = Databases(capture)
    DB.create_annot_table(dfannot)
    DB.create_doc_table(dfannot)


def get_annot_bedlocation(capture, pipelinedir=None):
    """Find annotated BEDfile in pipeline home. Return a filelocation."""
    bed_annot = os.path.join(pipelinedir, 'captures',
                             '{}_target.annotated'.format(capture))
    if not os.path.isfile(bed_annot):
        raise IOError('{} does not exist.'.format(bed_annot))
    return bed_annot


def add_docfile(docfile, capture, serie, sample):
    D = Databases(capture)
    data = D.parse_docfile(docfile)
    D.add_data_to_db(sample, serie, data)


def drop_badsamples(df, badsamples):
    drop = [i for i in badsamples if i in df.index.droplevel(0)]
    dfclean = df.drop(drop, level=1)
    return dfclean


def get_poscondata(df, posconsamples, keepseries=False):
    dfall = df.copy()
    if not keepseries:
        dfall.index = dfall.index.droplevel(0)
    try:
        df_poscon = dfall.loc[posconsamples]
    except KeyError:
        return pd.DataFrame()

    return df_poscon


def collect_archive(capture, correctmales=False):
    """Call Database class. Return dataframe"""
    archive = Databases(capture).get_archive()
    if correctmales:
        archive = correct_males(archive)
    return archive


def clean_archive(archive, capture, keepseries=False, getposcondf=False):
    QD = Databases(capture)
    pcsamples = QD.get_positive_controls()
    badsamples = QD.get_bad_samples()
    if badsamples:
        archive = drop_badsamples(archive, badsamples)
    if pcsamples:
        df_poscons = get_poscondata(archive, pcsamples, keepseries=keepseries)
        archive.drop(pcsamples, level=1, inplace=True)
    else:
        df_poscons = pd.DataFrame()
    if not keepseries:
        archive.index = archive.index.droplevel(0)
    if getposcondf:
        return (archive, df_poscons)
    else:
        return archive


def seperate_data(df, new, sample=None, keepseries=False):
    """Seperate dataframe in archive and to be analyzed (new) data.
    Return 2 dataframes.
    """
    if sample is not None:
        df_new = df.loc[new]
        try:
            df_new = df_new.loc[sample]
        except KeyError:
            print('Requested sample not in serie.')
            raise
        else:
            df.drop([sample], level=1, inplace=True)
    else:
        df_new = df.loc[new]
        df.drop([new], level=0, inplace=True)
    return (df_new, df)


def get_intervals_for_gene(annot, gene):
    """Get intervals for gene from annotdf. Returns (list/index)"""
    return annot[annot['gen'].str.contains('^{}$'.format(gene))].index


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


def get_gene_list(genefile):
    with open(genefile, 'r') as f:
        return [line.strip() for line in f]


def analyze(capture, serie, docfile=None, sample=None,
            outdir=None, genelist=None, addonly=False):
    if docfile:
        add_docfile(docfile, capture, serie, sample)
        if addonly:
            sys.exit()

    df = collect_archive(capture, correctmales=True)

    df_new, df_archive = seperate_data(df, serie, sample)

    df_archive, df_poscons = clean_archive(df_archive, capture,
                                           getposcondf=True)

    archtargetmean, archtargetstd = get_target_info(df_archive.transpose())
    targetinfo = archtargetmean.join(archtargetstd)

    try:
        df_archive.index = df_archive.index.droplevel(0)
    except AttributeError:
        pass

    """For each sample/poscon that needs to be analyzed:
    1) Transpose df to get targets as index
    2) Join archive-df with sample-df on index.
    3) Get z-scores for the joined dataframe
    4) Concat all the scores into 1 dataframe with list comprehension
    """

    if not df_poscons.empty:
        samples = df_poscons.index
        zscores_poscons = pd.concat(
            [get_zscore_df(df_archive.transpose().join(
             df_poscons.transpose()[sample]))[sample] for sample in samples],
            axis=1)

    if sample:
        samples = sample
        zscores_sample = get_zscore_df(df_archive.transpose().join(
            df_new.transpose()))[sample]

    elif serie and not sample:
        samples = df_new.index
        zscores_sample = pd.concat(
            [get_zscore_df(df_archive.transpose().join(
             df_new.transpose()[sample]))[sample] for sample in samples],
            axis=1)
    zscore_archive = get_zscore_df(df_archive.transpose())

    if not df_poscons.empty:
        data = zscores_poscons.join(zscore_archive.join(zscores_sample))

    elif df_poscons.empty:
        data = zscore_archive.join(zscores_sample)

    QD = Databases(capture)
    poscons = QD.get_positive_controls_dict()
    badsamples = QD.get_bad_samples()
    annot = QD.get_annot()
    empiricalfragments = QD.get_regions_to_exclude()
    config = get_config_dict('{}/config.py'.format(
        os.path.dirname(os.path.abspath(__file__))))

    if not outdir:
        outdir = config['outputdir']
        newdir = create_dirs(None, capture, serie, outdir)
    elif outdir:
        newdir = create_dirs(outdir, capture, serie, outdir)

    with open('{}/archive.txt'.format(newdir), 'w') as f:
        [f.write('{}\n'.format(sample))
         for sample in df_archive.index]

    badregions = get_bad_regions(df_archive.transpose())
    badregions = badregions.join(annot)

    if not badregions.empty:
        badregions.index = remove_underscore_targets(badregions.index)
        badregions.to_csv('{}/excluded.txt'.format(newdir), sep='\t')
        badregions.index = underscore_targets(badregions.index)

    with open('{}/excluded.txt'.format(newdir), 'a') as f:
        [[f.write('{}\tEmpirisch bepaald\n'.format(x)) for x in i.split(', ')]
         for i in empiricalfragments]
        if not badregions.empty:
            f.write('\n{0:0.1f} % is niet callable'.format(
                percentage_callable(annot.reset_index(),
                                    badregions.reset_index()) * 100))
        if badregions.empty:
            f.write('Geen.')

    dfclean = df_new.transpose()
    [dfclean.drop(i, axis=1, inplace=True) for i in badsamples if i in samples]
    dfnewmean, dfnewstd = get_target_info(dfclean)
    dfarchmean, dfarchstd = get_target_info(df_archive.transpose())

    pdf = PdfPages('{}/{}.pdf'.format(newdir, serie))
    Plotter = SeriePlots(capture, serie, pdf)

    Plotter.plot_qc_serie(dfnewmean.join(dfnewstd),
                          dfarchmean.join(dfarchstd),
                          len(df_archive.index))
    pdf.close()
    print('{} QC done'.format(serie))

    for i, sample in enumerate(samples):
        df_new_norm = normalize_df(df_new.loc[sample])
        sampleplotdf = annot.join(archtargetmean.join(df_new_norm))
        SaPlotter = SamplePlots(sample, newdir, badsamples=badsamples)
        SaPlotter.sample_qc(sampleplotdf)
        if not df_poscons.empty:
            data = zscores_poscons.join(
                zscore_archive.join(zscores_sample[sample]))

        elif df_poscons.empty:
            data = zscore_archive.join(zscores_sample[sample])
        # data[(data < -3) | (data > 3)].to_csv('out.csv')
        calls = data[(data[sample] < -3) | (data[sample] > 3)]

        if len(calls.index) != 0:
            if isinstance(genelist, list):
                reportgenes = list(genelist)
            else:
                reportgenes = get_gene_list(genelist)
            calls = calls.transpose()
            calls_per_target = {target: '{}/{}'.format(len(calls[(calls[target] > 3) | (calls[target] < -3)][target]),
                                                       len(calls[target]))
                                for target in calls}

            calls = calls.transpose()
            calls = annot.join(calls, how='right')
            genes = calls['gen'].unique()
            calls = calls[[sample, 'gen']].join(pd.DataFrame.from_dict(
                calls_per_target, orient='index')).join(badregions[['Mean', 'Std']]).fillna('OK')
            calls.columns = ['Z-score', 'Gen', 'Freq', 'Mean', 'Std']
            calls.index = remove_underscore_targets(calls.index)
            calls.to_csv('{}/Calls/{}.txt'.format(newdir, sample), sep='\t')
            samplepdf = PdfPages('{}/Calls/{}.pdf'.format(newdir, sample))

            for gene in genes:
                if genelist:
                    if gene not in reportgenes:
                        continue
                intervalstoplot = get_intervals_for_gene(annot, gene)
                datatoplot = filter_data_by_intervals(data, intervalstoplot)
                datatoplot = sort_by_interval(datatoplot.copy())
                targetinfo_toplot = filter_data_by_intervals(targetinfo,
                                                             intervalstoplot)
                targetinfo_toplot = sort_by_interval(targetinfo_toplot.copy())
                SaPlotter.plot_cnv_calls(datatoplot, gene, samplepdf,
                                         targetinfo_toplot, poscons)

            samplepdf.close()
        print('{} ({} of {}) done'.format(sample, i + 1, len(samples)))
