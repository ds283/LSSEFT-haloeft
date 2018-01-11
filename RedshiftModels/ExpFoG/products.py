import numpy as np

import sqlite3


# container class for theory products associated with an exponential (Gaussian) fingers-of-God model
class database(object):

    def __init__(self, my_config, k_sample, col_name):

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
                       'b1_1_bGamma3'], theory_db, my_config, col_name)


    def __import(self, tables, db, my_config, col_name):

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

            # for each power spectrum table, read in its mu0, mu2, mu4, mu6 and mu8 values
            for tag in tables:

                self.payload[tag] = self.__import_Pk(conn, data, tag, col_name)

            # need f to compute mu^6 counterterm, so read its value
            self.__import_f(conn, data)


    def __import_Pk(self, conn, data, tag, col_name):

        mu0 = self.__import_mu(conn, tag, data, 0, col_name)
        mu2 = self.__import_mu(conn, tag, data, 2, col_name)
        mu4 = self.__import_mu(conn, tag, data, 4, col_name)
        mu6 = self.__import_mu(conn, tag, data, 6, col_name)
        mu8 = self.__import_mu(conn, tag, data, 8, col_name)

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __import_mu(self, conn, tag, data, n, col_name):

        # obtain a database cursor
        cursor = conn.cursor()

        # construct the relevant table name from knowing its tag, and whether we want the n=0, 2, 4, 6, or 8th power
        table_name = '{tag}_mu{power}'.format(tag=tag, power=n)

        # execute SQL query
        cursor.execute(
            ("SELECT k_config.k AS k, sample." + col_name + " AS Pn FROM (SELECT * FROM " + table_name + " "
             "WHERE mid=:model AND growth_params=:growth AND loop_params=:loop "
             "AND zid=:zid AND init_Pk_id=:init_Pk AND final_Pk_id=:final_Pk AND IR_id=:IR_cutoff "
             "AND UV_id=:UV_cutoff) AS sample "
             "INNER JOIN k_config ON sample.kid = k_config.id ORDER BY k;"),
            data)

        ks = []
        Pns = []

        # read results from cursor
        for row in cursor:

            k, Pn = row
            ks.append(k)
            Pns.append(Pn)

        return np.asarray(Pns)


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
