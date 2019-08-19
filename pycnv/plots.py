#!/usr/bin/env python
import os
import math

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style('darkgrid')

class SeriePlots(object):
    """Create QC plot per serie."""

    def __init__(self, capture, serie, pdf):
        """Set capture, serie and PDF."""
        self.capture = capture
        self.serie = serie
        self.pdf = pdf

    def plot_qc_serie(self, df_new, df_arch, nrarchive=None):
        """Plot mean serie vs mean archive and stddev serie and archive."""
        fig = plt.figure(figsize=(12, 9))
        ax = plt.subplot(211)
        plt.title('{} QC {}'.format(self.serie, self.capture),
                  size=20)

        ax.scatter(df_arch['Mean'], df_new['Mean'], label=self.serie)
        ax.set_xlim([0, 3.5])
        ax.set_ylim([0, 3.5])
        ax.set_title('Gemiddelde genormaliseerde coverage per target',
                     size=20)
        if nrarchive is not None:
            ax.set_xlabel('Archief: {} samples'.format(nrarchive),
                          size=15)
        ax.set_ylabel('{}'.format(self.serie), size=15)
        ax.legend()

        ax2 = plt.subplot(223)
        ax2.hist(df_arch['Std'].dropna().values, 150)
        ax2.set_xlim([0, 0.5])
        ax2.set_title('Standaarddeviatie archief')

        ax3 = plt.subplot(224)
        ax3.hist(df_new['Std'].dropna().values, 150)
        ax3.set_xlim([0, 0.5])
        ax3.set_title('Standaarddeviatie nieuwe serie')

        fig.tight_layout()
        self.pdf.savefig()
        plt.close()


class SamplePlots(object):

    def __init__(self, sample, outdir=os.getcwd(), badsamples=[]):

        self.sample = sample
        self.outdir = outdir
        if self.sample in badsamples:
            self.badsample = True
        if self.sample not in badsamples:
            self.badsample = False

    def sample_qc(self, df, serie):
        axmax = list()
        axmin = list()

        genes = df.gen.unique()

        fig = plt.figure(figsize=(12, 9))
        ax = plt.subplot(111)

        colorlist = sns.color_palette("husl", len(genes))

        for i, g in enumerate(genes):
            intervalstoplot = df[df['gen'] == g]
            if intervalstoplot.empty:
                continue
            x = list(intervalstoplot['Mean'].values)
            y = list(intervalstoplot[self.sample].values)
            axmax.append(max(x+y))
            axmin.append(min(x+y))

            try:
                ax.scatter(x, y, label=g, color=colorlist[i])
            except ValueError as e:
                ax.scatter(x, y, label=g, color='grey')
                print(e, i)

        if len(genes) > 100:
            ax.legend(ncol=2, loc='center left', fontsize=2)
        else:
            ax.legend(ncol=2, loc='center left')
        axmax = math.ceil(max(axmax))
        axmin = math.floor(min(axmin))

        if axmax < 2:
            axmax = 2
        if axmin > 0:
            axmin = 0

        ax.set_xlim([axmin, axmax])
        ax.set_ylim([axmin, axmax])
        a = [_ for _ in range(axmin, axmax+1, 1)]
        b = [_ for _ in range(axmin, axmax+1, 1)]
        ax.plot(a, b)
        if self.badsample:
            plt.title('{} nonarchive'.format(self.sample), size=15)
        elif not self.badsample:
            plt.title(self.sample, size=15)
        fig.tight_layout()
        plt.savefig('{}/QC/{}.png'.format(self.outdir, self.sample),
                    dpi=80)
        plt.close()

    def plot_cnv_calls(self, data, gene, pdf, targetinfo, serie, poscons=None):
        if poscons is None:
            poscons = dict()
        fig = plt.figure(figsize=(12, 9))
        ax = plt.subplot2grid((8, 1), (0, 0), rowspan=4)
        ax2 = plt.subplot2grid((8, 1), (4, 0), rowspan=2, sharex=ax)
        ax3 = plt.subplot2grid((8, 1), (6, 0), rowspan=2, sharex=ax)

        for pat in data:
            y = data[pat].values
            x = np.arange(1, len(y) + 1, 1.0)

            if pat == self.sample:
                if self.badsample:
                    ax.plot(x, y, 'ro', markersize=8,
                            label='{}\nnonarchive'.format(pat))
                elif not self.badsample:
                    ax.plot(x, y, 'ro', markersize=8, label=pat)
            elif pat in serie:
                ax.plot(x, y, 'bo', markersize=4)
            # Pos control AND correct gene => label = PosConDnr + Gene
            elif pat in poscons and gene in poscons[pat]:
                ax.plot(x, y, '--', markersize=4, label=('PosCon'))

            elif pat in poscons and gene not in poscons[pat]:
                continue
            # Archive = black label
            else:
                ax.plot(x, y, 'ko', markersize=4)

        if max(data.max().sort_values()) > 10:
            maxax = max(data.max().sort_values())
            maxax += 0.5
        else:
            maxax = 10
        if min(data.min().sort_values()) < -10:
            minax = min(data.min().sort_values())
            minax -= 0.5
        else:
            minax = -10

        ax.set_ylim([minax, maxax])
        ax.set_ylabel('Z-score', size=15)
        ax.set_title(gene, size=20)
        ax.axhline(y=3)
        ax.axhline(y=-3)

        z = targetinfo['Std'].values
        x = range(1, len(z) + 1)
        ax2.plot(x, z, 'ko', markersize=6)
        ax2.axhline(y=0.15, c='k')
        ax2.set_ylim([0, 0.5])
        ax2.set_ylabel('Stdev. genorm. coverage', size=8)
        ax2.set_title('Variatie')

        z = targetinfo['Mean'].values
        x = range(1, len(z) + 1)
        ax3.plot(x, z, 'bo', markersize=6)
        ax3.axhline(y=0.2, c='b')
        ax3.set_ylim([0, 2])
        ax3.set_ylabel('Genorm. coverage', size=8)
        ax3.set_title('Gemiddelde per regio')

        xticks = np.arange(0, len(data) + 2, 1.0)

        ax.xaxis.set_ticks(xticks)
        ax2.xaxis.set_ticks(xticks)
        ax3.xaxis.set_ticks(xticks)

        ax.legend(fontsize=8, ncol=1, loc='upper right',
                  bbox_to_anchor=(1.05, 1))

        fig.tight_layout(pad=4)
        pdf.savefig()

        plt.close()
