import numpy as np
from collections import OrderedDict

import sqlite3


# container class for EFT theory products
class database(object):


    def __init__(self, my_config, k_sample):

        self.k_sample = k_sample


        # READ AND CACHE THEORY DATA PRODUCTS

        theory_db = my_config["HaloEFT", "theory_db"]

        self.payload = {}

        # map from my table labels to Lucia's table names
        tables = OrderedDict([('nobias', "One"),
                              ('b1_1', "b1o1"),
                              ('b1_2', "b1o2"),
                              ('b1_3', "b1o3"),
                              ('b2_2', "b2o2"),
                              ('b2_3', "b2o3"),
                              ('b3', "b3"),
                              ('bG2_2', "bg2o2"),
                              ('bG2_3', "bg2o3"),
                              ('bdG2', "b1g2"),
                              ('bGamma3', "bGam3"),
                              ('b1_1_b1_1', "b1o1b1o1"),
                              ('b1_2_b1_2', "b1o2b1o2"),
                              ('b1_1_b1_2', "b1o2b1o2"),
                              ('b1_1_b1_3', "b1o1b1o3"),
                              ('b1_1_b2_2', "b1o1b2o2"),
                              ('b1_1_b2_3', "b1o1b2o3"),
                              ('b1_2_b2_2', "b1o2b2o2"),
                              ('b2_2_b2_2', "b2o2b2o2"),
                              ('b1_1_b3', "b1o1b3"),
                              ('b1_1_bG2_2', "b1o1bg2o2"),
                              ('b1_1_bG2_3', "b1o1bg2o3"),
                              ('b1_2_bG2_2', "b1o2bg2o2"),
                              ('bG2_2_bG2_2', "bg2o2bg2o2"),
                              ('b2_2_bG2_2', "b2o2bg2o2"),
                              ('b1_1_bdG2', "b1o1b1g2"),
                              ('b1_1_bGamma3', 'b1o1bGam3')])

        self.__import(tables, ['c0', 'c2', 'c4', 'c6'], theory_db, my_config)


    def __import(self, tables, counterterms, db, my_config):

        # open SQLite3 connexion to database
        with sqlite3.connect(db) as conn:

            # for each power spectrum table, read in its P0, P2, and P4 values
            for tag in tables:

                self.payload[tag] = self.__import_Pk(conn, tables[tag])

            # for each counterterm, read in its values likewise
            for tag in counterterms:

                self.payload[tag] = self.__import_counterterm(conn, tag)

            # need f to compute mu^6 counterterm, so read its value
            self.__import_f()

        # finally, construct stochastic counterterms
        ks = self.k_sample.conv_ks
        ksq = ks*ks

        d1_P0 = np.power(ks, 0)
        d1_P2 = 0*ks
        d1_P4 = 0*ks

        d2_P0 = ksq
        d2_P2 = 0*ks
        d2_P4 = 0*ks

        d3_P0 = ksq/3.0
        d3_P2 = 2.0*ksq/3.0
        d3_P4 = 0*ks

        self.payload['d1'] = np.array([d1_P0, d1_P2, d1_P4])
        self.payload['d2'] = np.array([d2_P0, d2_P2, d2_P4])
        self.payload['d3'] = np.array([d3_P0, d3_P2, d3_P4])


    def __import_Pk(self, conn, table_id):

        P0 = self.__import_P_ell(conn, table_id, 0)
        P2 = self.__import_P_ell(conn, table_id, 2)
        P4 = self.__import_P_ell(conn, table_id, 4)

        Pk_group = np.array([P0, P2, P4])

        return Pk_group


    def __import_P_ell(self, conn, table_id, ell):

        # obtain a database cursor
        cursor = conn.cursor()

        # construct the relevant table name from knowing its tag, and whether we want the ell=0, 2, or 4 mode
        table_name = 'P{ell}'.format(ell=ell)
        column_name = 'T{ell}{id}'.format(ell=ell, id=table_id)

        # execute SQL query
        cursor.execute("SELECT k AS k, " + column_name + " AS Pell FROM " + table_name + " ORDER BY k;")

        ks = []
        Pells = []

        # read results from cursor
        for row in cursor:

            k, Pell = row
            ks.append(k)
            Pells.append(Pell)

        return np.asarray(Pells)


    def __import_counterterm(self, conn, tag):

        # obtain a database cursor
        cursor = conn.cursor()

        # construct the relevant table name from knowings its tag
        table_name = 'counterterms_{tag}'.format(tag=tag)

        # execute SQL query
        cursor.execute("SELECT k AS k, P0 AS P0, P2 AS P2, P4 AS P4 FROM " + table_name + " ORDER BY k;")

        ks = []
        P0s = []
        P2s = []
        P4s = []

        # read results from cursor
        for row in cursor:

            k, P0, P2, P4 = row
            ks.append(k)
            P0s.append(P0)
            P2s.append(P2)
            P4s.append(P4)

        return np.array([ np.asarray(P0s), np.asarray(P2s), np.asarray(P4s) ])


    def __import_f(self):

        self.f = 0.7039874052239977
