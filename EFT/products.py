import numpy as  np

import sqlite3


# container class for EFT theory products
class database(object):

    def __init__(self, my_config, k_sample):

        self.k_sample = k_sample


        # READ AND CACHE THEORY DATA PRODUCTS

        theory_db = my_config["HaloEFT", "theory_db"]

        self.payload = {}

        self.__import(['nobias',
                       'b1_1',
                       'b1_2',
                       'b1_3',
                       'b2_2',
                       'b2_3',
                       'b3',
                       'bG2_2',
                       'bG2_3',
                       'bdG2',
                       'bGamma3',
                       'b1_1_b1_1',
                       'b1_2_b1_2',
                       'b1_1_b1_2',
                       'b1_1_b1_3',
                       'b1_1_b2_2',
                       'b1_1_b2_3',
                       'b1_2_b2_2',
                       'b2_2_b2_2',
                       'b1_1_b3',
                       'b1_1_bG2_2',
                       'b1_1_bG2_3',
                       'b1_2_bG2_2',
                       'bG2_2_bG2_2',
                       'b2_2_bG2_2',
                       'b1_1_bdG2',
                       'b1_1_bGamma3'], ['c0', 'c2', 'c4', 'c6'], theory_db, my_config)


    def __import(self, tables, counterterms, db, my_config):

        # extract identifiers needed to fix which model we read from the database -- there are a lot of these,
        # needed to fix the final redshift, parameters used during integration, IR and UV cutoffs
        # and IR resummation scale
        model = my_config["HaloEFT", "model"]
        growth_params = my_config["HaloEFT", "growth_params"]
        loop_params = my_config["HaloEFT", "loop_params"]
        XY_params = my_config["HaloEFT", "XY_params"]
        zid = my_config["HaloEFT", "zid"]
        init_Pk = my_config["HaloEFT", "init_Pk"]
        final_Pk = my_config["HaloEFT", "final_Pk"]
        IR_cutoff = my_config["HaloEFT", "IR_cutoff"]
        UV_cutoff = my_config["HaloEFT", "UV_cutoff"]
        IR_resum = my_config["HaloEFT", "IR_resum"]

        # bundle identifiers together into a dictionary for easy use
        data = {'model': model, 'growth': growth_params, 'loop': loop_params, 'XY': XY_params,
                'zid': zid, 'init_Pk': init_Pk, 'final_Pk': final_Pk,
                'IR_cutoff': IR_cutoff, 'UV_cutoff': UV_cutoff, 'IR_resum': IR_resum}

        # open SQLite3 connexion to database
        with sqlite3.connect(db) as conn:

            # for each power spectrum table, read in its P0, P2, and P4 values
            for tag in tables:

                self.payload[tag] = self.__import_Pk(conn, data, tag)

            # for each counterterm, read in its values likewise
            for tag in counterterms:

                self.payload[tag] = self.__import_counterterm(conn, tag, data)

            # need f to compute mu^6 counterterm, so read its value
            self.__import_f(conn, data)

        # finally, construct stochastic counterterms
        ks = self.k_sample.WiggleZ_conv_ks
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


    def __import_Pk(self, conn, data, tag):

        P0 = self.__import_P_ell(conn, tag, data, 0)
        P2 = self.__import_P_ell(conn, tag, data, 2)
        P4 = self.__import_P_ell(conn, tag, data, 4)

        Pk_group = np.array([P0, P2, P4])

        return Pk_group


    def __import_P_ell(self, conn, tag, data, ell):

        # obtain a database cursor
        cursor = conn.cursor()

        # construct the relevant table name from knowing its tag, and whether we want the ell=0, 2, or 4 mode
        table_name = '{tag}_P{ell}'.format(tag=tag, ell=ell)

        # execute SQL query
        cursor.execute(
            ("SELECT k_config.k AS k, sample.P1loopSPT_resum AS Pell FROM (SELECT * FROM " + table_name + " "
             "WHERE mid=:model AND growth_params=:growth AND loop_params=:loop AND XY_params=:XY "
             "AND zid=:zid AND init_Pk_id=:init_Pk AND final_Pk_id=:final_Pk AND IR_cutoff_id=:IR_cutoff "
             "AND UV_cutoff_id=:UV_cutoff AND IR_resum_id=:IR_resum) AS sample "
             "INNER JOIN k_config ON sample.kid = k_config.id ORDER BY k;"),
            data)

        ks = []
        Pells = []

        # read results from cursor
        for row in cursor:

            k, Pell = row
            ks.append(k)
            Pells.append(Pell)

        return np.asarray(Pells)


    def __import_counterterm(self, conn, tag, data):

        # obtain a database cursor
        cursor = conn.cursor()

        # construct the relevant table name from knowings its tag
        table_name = 'counterterms_{tag}'.format(tag=tag)

        # execute SQL query
        cursor.execute(
            ("SELECT k_config.k AS k, sample.P0_k2_resum AS P0, sample.P2_k2_resum AS P2, sample.P4_k2_resum AS P4 "
             "FROM (SELECT * FROM " + table_name + " WHERE mid=:model AND growth_params=:growth AND XY_params=:XY "
             "AND zid=:zid AND init_Pk_id=:init_Pk AND final_Pk_id=:final_Pk AND IR_cutoff_id=:IR_cutoff "
             "AND UV_cutoff_id=:UV_cutoff AND IR_resum_id=:IR_resum) AS sample "
             "INNER JOIN k_config ON sample.kid = k_config.id ORDER BY k;"),
            data)

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


    def __import_f(self, conn, data):

        # obtain database cursor
        cursor = conn.cursor()

        # execute SQL query
        cursor.execute("SELECT f_linear FROM f_factors WHERE zid=:zid;", data)

        # read value
        flist = cursor.fetchone()

        if flist is None:
            raise LookupError

        self.f = flist[0]
