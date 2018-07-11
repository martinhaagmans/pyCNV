#!/usr/bin/env python
import sys
import os
import csv
import sqlite3
import json
import pandas as pd


class Databases:

    def __init__(self, capture, configfile=None, maxcol=2000):
        if configfile is None:
            configfile = '{}/config.py'.format(
                os.path.dirname(os.path.abspath(__file__)))
        self.capture = capture

        config = {}
        exec(open(configfile).read(), config)

        self.dir = config['dbdir']
        self.dbgeneral = config['dbgeneral']
        self.poscontable = config['poscontable']
        self.badsampletable = config['badsampletable']
        self.badregiontable = config['badregiontable']
        self.pipelinedir = config['pipelinedir']
        self.conn = sqlite3.connect('{}/{}'.format(self.dir, self.dbgeneral))
        self.capconn = sqlite3.connect('{}/{}.sqlite'
                                       .format(self.dir, self.capture))

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
samples
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

    @staticmethod
    def parse_docfile(docfile):
        target_coverage = list()
        with open(docfile) as f:
            fin = csv.reader(f, delimiter='\t')
            _ = next(fin)
            for line in fin:
                target, total, mean, *_ = line
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
        except sqlite3.IntegrityError as e:
            print(e)
        else:
            self.capconn.commit()

    def create_badregiontable(self):
        sql = """CREATE TABLE IF NOT EXISTS {}
        (capture text NOT NULL,
        target text NOT NULL,
        gen text NOT NULL,
        PRIMARY KEY(capture, target))
        """.format(self.badregiontable)
        c = self.conn.cursor()
        c.execute(sql)
        self.conn.commit()

    def get_regions_to_exclude(self):
        c = self.conn.cursor()
        sql = "SELECT DISTINCT fragments FROM {} WHERE (capture='{}')".format(
            self.badregiontable, self.capture)
        c.execute(sql)
        fragments = [val for tup in c.fetchall() for val in tup]
        return fragments

    def create_annot_table(self, dfannot):
        """Create a table from the annotation dataframe."""

        checksql = """SELECT * FROM sqlite_master
        WHERE name ='{}annot' and type='table'
        """.format(self.capture)

        c = self.capconn.cursor()
        c.execute(checksql)

        if c.fetchall():
            print('Annot DB exists')

        else:
            sql = """CREATE TABLE IF NOT EXISTS {}annot
            (target text NOT NULL,
            gen text NOT NULL,
            PRIMARY KEY(target,gen))
            """.format(self.capture.lower())
            c.execute(sql)
            [c.execute("INSERT INTO {}annot VALUES ('{}', '{}')".format(
             self.capture, i, dfannot.loc[i]['Gen']))
             for i in dfannot.index.unique()]
            self.capconn.commit()

    def create_numbered_table(self, targets, tablenumber, c):
        createsql = """CREATE TABLE IF NOT EXISTS {}{}
        (sample text NOT NULL,
        serie text NOT NULL,
        PRIMARY KEY(sample,serie))
        """.format(self.capture, tablenumber)

        c.execute(createsql)
        [c.execute("ALTER TABLE {}{} ADD COLUMN '{}' REAL".format(
         self.capture, tablenumber, i)) for i in targets]

    def create_doc_table(self, dfannot):
        """Create a table with 1 column per capture-target."""
        sql = '''CREATE TABLE IF NOT EXISTS {}
        (SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        DATA TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE))
        '''.format(self.capture.lower())
        c = self.capconn.cursor()
        c.execute(sql)
        self.capconn.commit()

    def add_poscontrols(self, inputfile):
        """Add poscontroles from file to table"""

        sql = """CREATE TABLE IF NOT EXISTS {}
        (capture text NOT NULL,
        gene text NOT NULL,
        sample text NOT NULL,
        soort text NOT NULL,
        PRIMARY KEY(capture, gene,sample))
        """.format(self.poscontable)

        c = self.conn.cursor()
        c.execute(sql)

        with open(inputfile, 'r') as f:
            for line in f.read().splitlines():
                if not line.strip():
                    continue
                cap, gen, sample, soort = (line.split())
                sql = ("INSERT INTO {} VALUES ('{}', '{}', '{}', '{}')"
                       .format(self.poscontable, cap, gen,
                               sample, soort))
                try:
                    c.execute(sql)
                except sqlite3.IntegrityError as e:
                    print(e)
        self.conn.commit()

    def add_badsamples(self, inputfile):
        """Add badsamples from file to table"""

        sql = """CREATE TABLE IF NOT EXISTS {}
        (serie text NOT NULL,
        sample text NOT NULL,
        code text,
        PRIMARY KEY(serie,sample))
        """.format(self.badsampletable)

        c = self.conn.cursor()
        c.execute(sql)

        with open(inputfile, 'r') as f:
            for line in f.read().splitlines():
                if not line.strip():
                    continue
                serie, sample, code = line.split()[:3]
                if 'Serie' not in serie:
                    sql = ("INSERT INTO {} VALUES ('Serie{}', '{}', "
                           "'{}')".format(self.badsampletable, serie,
                                          sample, code))
                else:
                    sql = ("INSERT INTO {} VALUES ('{}', '{}', '{}')"
                           .format(self.badsampletable, serie,
                                   sample, code))
                try:
                    c.execute(sql)
                except sqlite3.IntegrityError as e:
                    print(e)
        self.conn.commit()

    def get_bad_samples(self):
        """Get badsample ID's from table and return a list."""

        checksql = """SELECT * FROM sqlite_master
        WHERE name ='{}' and type='table'
        """.format(self.badsampletable)

        c = self.conn.cursor()
        c.execute(checksql)

        if not c.fetchall():
            print('No badsamples table detected')
            return []
        else:
            sql = 'SELECT DISTINCT sample FROM {} '.format(self.badsampletable)
            c.execute(sql)
            badsamples = [val for tup in c.fetchall() for val in tup]
            return badsamples

    def get_positive_controls(self):
        """Get positive control info from table and return a list."""

        checksql = """SELECT * FROM sqlite_master
        WHERE name ='{}' and type='table'
        """.format(self.poscontable)

        c = self.conn.cursor()
        c.execute(checksql)

        if not c.fetchall():
            print('No positive control samples table detected')
            return []
        else:
            sql = "SELECT DISTINCT sample FROM {} WHERE (capture='{}')".format(
                self.poscontable, self.capture)
            c.execute(sql)
            poscons = [val for tup in c.fetchall() for val in tup]
            return poscons

    def get_positive_controls_dict(self):
        """Get positive control info from table and return a dict."""

        checksql = """SELECT * FROM sqlite_master
        WHERE name ='{}' and type='table'
        """.format(self.poscontable)

        c = self.conn.cursor()
        c.execute(checksql)

        if not c.fetchall():
            print('No positive control samples table detected')
            return []
        else:
            sql = ("SELECT DISTINCT sample, gene FROM {} WHERE (capture='{}')"
                   .format(self.poscontable, self.capture))
            c.execute(sql)
            poscons = {sample: gene for (sample, gene) in c.fetchall()}
            return poscons

    def get_annot(self):
        """Get annotation(target-gene) from table and return a df."""
        sql = "SELECT * FROM {}annot".format(self.capture)

        try:
            annot = pd.read_sql(sql, self.capconn)
        except pd.io.sql.DatabaseError as e:
            print(e)
            sys.exit()
        else:
            annot.set_index('target', inplace=True, drop=True)
            return annot

    def correct_posconserie(self):
        """Correct serie for positive control into PosCon."""

        checksql = """SELECT * FROM sqlite_master
        WHERE name ='{}' and type='table'
        """.format(self.poscontable)

        c = self.conn.cursor()
        c.execute(checksql)

        if not c.fetchall():
            print('No positive control table detected')
            return

        else:
            sql = ("SELECT DISTINCT sample, capture FROM {}"
                   .format(self.poscontable))
            c.execute(sql)
            poscons = {capture: sample for (sample, capture) in c.fetchall()}

        for cap in poscons.keys():
            db = '{}/{}.sqlite'.format(self.dir, cap)
            conn = sqlite3.connect(db)
            cc = conn.cursor()
            cc.execute(checksql)

            if not c.fetchall():
                print('No {} table detected'.format(cap))
                continue
            else:
                [cc.execute("UPDATE {} SET serie='PosCon' "
                            "WHERE (sample='{}' AND serie!='PosCon')"
                            .format(cap, i)) for i in poscons[cap]]
            conn.close()

        def remove(sample=None, serie=None):
            if sample and serie is None:
                sql = '''DELETE FROM {}
                WHERE SAMPLE='{}'
                '''.format(self.capture.lower(), sample)

            elif serie and sample is None:
                sql = '''DELETE FROM {}
                WHERE SERIE='{}'
                '''.format(self.capture.lower(), serie)

            elif sample and serie:
                sql = '''DELETE FROM {}
                WHERE (SAMPLE='{}' AND SERIE='{}')
                '''.format(self.capture.lower(), sample, serie)

            elif sample is None and serie is None:
                raise ValueError('Geen sample en/of serie opgegeven')

            c = self.capconn.cursor()
            c.execute(sql)
            self.capconn.commit()


class UpdateDatabase(object):

    def __init__(self, new, old, configfile=None):
        if configfile is None:
            configfile = '{}/config.py'.format(
                os.path.dirname(os.path.abspath(__file__)))

        config = {}
        exec(open(configfile).read(), config)

        self.poscontable = config['poscontable']
        self.badsampletable = config['badsampletable']
        self.badregiontable = config['badregiontable']
        self.conn = sqlite3.connect('{}/{}'.format(config['dbdir'],
                                                   config['dbgeneral']))
        self.c = self.conn.cursor()

        self.new = new
        self.old = old

    def poscons(self):
        sql = '''INSERT INTO {tn} (capture, gene, sample, soort)
                    SELECT "{newcap}" , gene, sample, soort
                        FROM {tn} WHERE capture="{oldcap}"
        '''.format(tn=self.poscontable, newcap=self.new, oldcap=self.old)
        self.c.execute(sql)
        self.conn.commit()

    def badregions(self):
        sql = '''INSERT INTO {tn} (capture, fragments)
                    SELECT "{newcap}", fragments
                        FROM {tn} WHERE capture="{oldcap}"
        '''.format(tn=self.badregiontable, newcap=self.new, oldcap=self.old)
        self.c.execute(sql)
        self.conn.commit()

    def all(self):

        try:
            self.poscons()
        except sqlite3.IntegrityError as e:
            print(e)

        try:
            self.badregions()
        except sqlite3.IntegrityError as e:
            print(e)

        self.conn.close()
