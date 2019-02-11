#!/usr/bin/env python

import os
import sys
import csv
import json
import sqlite3

import pandas as pd


class Databases:

    def __init__(self, capture, configfile=None, maxcol=2000):

        scriptdir = os.path.dirname(os.path.abspath(__file__))

        if configfile is None:
            configfile = os.path.join(scriptdir, 'config.py')
        self.capture = capture
        config = dict()
        exec(open(configfile).read(), config)
        db_dir = config['dbdir']
        db_general = os.path.join(db_dir, config['dbgeneral'])
        db_capture = os.path.join(db_dir, '{}.sqlite'.format(capture))
        self.table_poscons = config['poscontable']
        self.table_badsamples = config['badsampletable']
        self.table_badregion = config['badregiontable']
        self.conn = sqlite3.connect(db_general)
        self.capconn = sqlite3.connect(db_capture)

    def create_bad_region_table(self):
        """Create a table to hold empirical exluded regions."""
        sql = """CREATE TABLE IF NOT EXISTS {}
        (capture text NOT NULL,
        target text NOT NULL,
        gen text NOT NULL,
        PRIMARY KEY(capture, target))
        """.format(self.table_badregion)
        c = self.conn.cursor()
        c.execute(sql)
        self.conn.commit()

    def create_bad_sample_table(self):
        """Create a table to hold samples to be exluded from archive."""
        sql = """CREATE TABLE IF NOT EXISTS {}
        (serie text NOT NULL,
        sample text NOT NULL,
        code text,
        PRIMARY KEY(serie,sample))
        """.format(self.table_badsamples)
        c = self.conn.cursor()
        c.execute(sql)
        self.conn.commit()

    def create_annot_table(self, dfannot):
        """Create a table from the annotation dataframe."""
        sql = """CREATE TABLE IF NOT EXISTS {}annot
        (target text NOT NULL,
        gen text NOT NULL,
        PRIMARY KEY(target, gen))
        """.format(self.capture.lower())
        c = self.capconn.cursor()
        c.execute(sql)
        [c.execute("INSERT INTO {}annot VALUES ('{}', '{}')".format(
            self.capture, i, dfannot.loc[i]['Gen']))
            for i in dfannot.index.unique()]
        self.capconn.commit()
        return

    def create_doc_table(self):
        """Create a table to hold DoC data."""
        sql = '''CREATE TABLE IF NOT EXISTS {}
        (SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        DATA TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE))
        '''.format(self.capture.lower())
        c = self.capconn.cursor()
        c.execute(sql)
        self.capconn.commit()
        return

    def add_poscontrols(self, inputfile):
        """Add poscontroles from file to table"""
        table = self.table_poscons
        c = self.conn.cursor()

        with open(inputfile, 'r') as f:
            for line in f:
                if not line:
                    continue
                cap, gen, sample, soort = line.split()
                sql = """INSERT INTO {}
                VALUES ('{}', '{}', '{}', '{}')
                """.format(table, cap, gen, sample, soort)
                try:
                    c.execute(sql)
                except sqlite3.IntegrityError:
                    pass
        self.conn.commit()

    def add_badsamples(self, inputfile):
        """Add badsamples from file to table"""

        table = self.table_badsamples
        c = self.conn.cursor()

        with open(inputfile, 'r') as f:
            for line in f.read().splitlines():
                if not line.strip():
                    continue
                serie, sample, code = line.split()[:3]
                sql = """INSERT INTO {}
                VALUES ('{}', '{}', '{}')
                """.format(table, serie, sample, code)
                try:
                    c.execute(sql)
                except sqlite3.IntegrityError:
                    pass
        self.conn.commit()

    @staticmethod
    def parse_docfile(docfile):
        target_coverage = list()
        with open(docfile) as f:
            fin = csv.reader(f, delimiter='\t')
            _ = next(fin)
            for line in fin:
                target, _total, mean, *_ = line
                target = target.replace(':', '_')
                target = target.replace('-', '_')
                target_coverage.append((target, mean))
        return target_coverage

    def add_data_to_db(self, sample, serie, data):
        sql = """INSERT INTO {}
        (SAMPLE, SERIE, DATA)
        VALUES ('{}', '{}', '{}')
        """.format(self.capture.lower(), sample, serie, json.dumps(data))
        c = self.capconn.cursor()
        try:
            c.execute(sql)
        except sqlite3.IntegrityError:
            pass
        else:
            self.capconn.commit()
        return

    def get_regions_to_exclude(self):
        c = self.conn.cursor()
        sql = """SELECT DISTINCT fragments
        FROM {}
        WHERE (capture='{}')
        """.format(self.table_badregion, self.capture)
        c.execute(sql)
        fragments = [val for tup in c.fetchall() for val in tup]
        return fragments

    def get_bad_samples(self):
        """Get badsample ID's from table and return a list."""

        c = self.conn.cursor()
        sql = """SELECT DISTINCT sample
        FROM {}""".format(self.table_badsamples)
        c.execute(sql)
        badsamples = [val for tup in c.fetchall() for val in tup]
        return badsamples

    def get_positive_controls_dict(self):
        """Get positive control info from table and return a dict."""
        c = self.conn.cursor()
        sql = """SELECT DISTINCT sample, gene
        FROM {}
        WHERE (capture='{}')
        """.format(self.table_poscons, self.capture)
        c.execute(sql)
        poscons = {sample: gene for (sample, gene) in c.fetchall()}
        return poscons

    def get_annot(self):
        """Get annotation(target-gene) from table and return a df."""
        sql = """SELECT target, gen
        FROM {}annot""".format(self.capture.lower())
        annot = pd.read_sql(sql, self.capconn)
        annot.set_index('target', inplace=True, drop=True)
        return annot

    def get_archive(self):
        """Get coveragedata from table and return df"""
        sql = "SELECT * FROM {}".format(self.capture.lower())
        c = self.capconn.cursor()
        c.execute(sql)
        d = dict()
        for sample, serie, data in c.fetchall():
            d[sample] = dict()

            d[sample]['data'] = dict()
            d[sample]['serie'] = serie

            for tup in json.loads(data):
                target, coverage = tup
                d[sample]['data'][target] = coverage

        dflist = list()

        for sample, data in d.items():
            dftmp = pd.Series(d[sample]['data'], name=sample, dtype=float)
            dftmp['serie'] = d[sample]['serie']
            dftmp['sample'] = sample
            dftmp = dftmp.transpose()
            dflist.append(dftmp)

        df = pd.concat(dflist, axis=1, sort=True)
        df = df.transpose().set_index(['serie', 'sample']).sort_index().sort_index(axis=1)
        return df

    def delete_serie(self, serie):
        """Delete serie from database"""
        sql = """DELETE FROM {}
        WHERE SERIE='{}'
        """.format(self.capture.lower(), serie)
        c = self.capconn.cursor()
        c.execute(sql)
        self.capconn.commit()
        return
