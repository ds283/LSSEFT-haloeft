import numpy as np

from astropy.io import ascii
from scipy import optimize

import sqlite3


class HaloEFT_core(object):

    def __init__(self, my_config, my_name):

        self.__mod_name = my_name

        self.nbin = 15
        self.nbinc = 25

        self.WiggleZ_mean_ks = np.linspace(0.01, 0.29, self.nbin)
        self.WiggleZ_conv_ks = np.linspace(0.01, 0.49, self.nbinc)

        self.labels = ['0.01', '0.03', '0.05', '0.07', '0.09',
                       '0.11', '0.13', '0.15', '0.17', '0.19',
                       '0.21', '0.23', '0.25', '0.27', '0.29']

        # READ AND CACHE THEORY DATA PRODUCTS

        theory_db = my_config["HaloEFT", "theory_db"]

        self.theory_payload = {}

        self.__import_theory(['nobias',
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

        fit_kmin = my_config["HaloEFT", "fit_kmin"]
        fit_kmax = my_config["HaloEFT", "fit_kmax"]

        mam = np.all([self.WiggleZ_mean_ks > fit_kmin, self.WiggleZ_mean_ks <= fit_kmax], axis=0)
        cam = np.all([self.WiggleZ_conv_ks > fit_kmin, self.WiggleZ_conv_ks <= fit_kmax], axis=0)
        self.mean_all_mask = np.concatenate((mam, mam, mam))
        self.conv_all_mask = np.concatenate((cam, cam, cam))

        ren_kmin = my_config["HaloEFT", "renormalize_kmin"]
        ren_kmax = my_config["HaloEFT", "renormalize_kmax"]

        mrm = np.all([self.WiggleZ_mean_ks > ren_kmin, self.WiggleZ_mean_ks <= ren_kmax], axis=0)
        crm = np.all([self.WiggleZ_conv_ks > ren_kmin, self.WiggleZ_conv_ks <= ren_kmax], axis=0)
        self.__mean_ren_mask = np.concatenate((mrm, mrm, mrm))
        self.__conv_ren_mask = np.concatenate((crm, crm, crm))


        # READ AND CACHE WIGGLEZ DATA PRODUCTS

        h1_means = my_config["HaloEFT", "h1_means"]
        h1_cov = my_config["HaloEFT", "h1_matrix"]

        h3_means = my_config["HaloEFT", "h3_means"]
        h3_cov = my_config["HaloEFT", "h3_matrix"]

        h9_means = my_config["HaloEFT", "h9_means"]
        h9_cov = my_config["HaloEFT", "h9_matrix"]

        h11_means = my_config["HaloEFT", "h11_means"]
        h11_cov = my_config["HaloEFT", "h11_matrix"]

        h15_means = my_config["HaloEFT", "h15_means"]
        h15_cov = my_config["HaloEFT", "h15_matrix"]

        h22_means = my_config["HaloEFT", "h22_means"]
        h22_cov = my_config["HaloEFT", "h22_matrix"]

        self.data_raw_means = {}
        self.data_all_means = {}
        self.data_ren_means = {}
        self.data_all_covs = {}
        self.data_ren_covs = {}
        self.data_convs = {}
        self.data_regions = ['1h', '3h', '8h', '11h', '15h', '22h']

        self.__import_WiggleZ_data("1h", h1_means, h1_cov, my_config)
        self.__import_WiggleZ_data("3h", h3_means, h3_cov, my_config)
        self.__import_WiggleZ_data("8h", h9_means, h9_cov, my_config)
        self.__import_WiggleZ_data("11h", h11_means, h11_cov, my_config)
        self.__import_WiggleZ_data("15h", h15_means, h15_cov, my_config)
        self.__import_WiggleZ_data("22h", h22_means, h22_cov, my_config)


    def __import_WiggleZ_data(self, tag, means_file, matrix_file, my_config):

        realization = my_config["HaloEFT", "realization"]

        # number of bins in P_ell measurements and covariance matrix
        nbin = self.nbin

        # number of bins in convolution  matrix
        nbinc = self.nbinc

        # read in the WiggleZ data for this region as astropy.table.Table instances
        mean_table = ascii.read(means_file, Reader=ascii.NoHeader,
                                names=['ireal', 'k', 'P0', 'P0err', 'P2', 'P2err', 'P4', 'P4err'])

        cov_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=0, data_end=(3*nbin)**2,
                               names=['i', 'j', 'value'])

        conv_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=(3*nbin)**2+1, data_end=(3*nbin)**2+(3*nbinc)**2,
                                names=['i', 'j', 'value'])

        # want to extract the power spectrum for a given realization, so group mean_table by realization
        mean_by_ireal = mean_table.group_by('ireal')
        realization_mask = mean_by_ireal.groups.keys['ireal'] == realization
        mean_realization = mean_by_ireal.groups[realization_mask]

        # sort into ascending order of k
        mean_realization.sort('k')

        this_P0 = np.asarray(mean_realization['P0'])
        this_P2 = np.asarray(mean_realization['P2'])
        this_P4 = np.asarray(mean_realization['P4'])

        self.data_raw_means[tag] = np.concatenate((this_P0, this_P2, this_P4))
        self.data_all_means[tag] = np.concatenate((this_P0, this_P2, this_P4))[self.mean_all_mask]
        self.data_ren_means[tag] = np.concatenate((this_P0, this_P2, this_P4))[self.__mean_ren_mask]

        CMat = np.empty((3*nbin, 3*nbin))
        ConvMat = np.empty((3*nbinc, 3*nbinc))

        for row in cov_table:

            CMat[row['i']-1, row['j']-1] = row['value']

        for row in conv_table:

            ConvMat[row['i']-1, row['j']-1] = row['value']

        CMat_all = (CMat[self.mean_all_mask, :])[:, self.mean_all_mask]

        w, p = np.linalg.eig(CMat_all)
        if not np.all(w > 0):
            print 'using pseudo-inverse covariance matrix for "all" group in region {tag}'.format(tag=tag)
            self.data_all_covs[tag] = np.linalg.pinv(CMat_all)
        else:
            self.data_all_covs[tag] = np.linalg.inv(CMat_all)

        CMat_ren = (CMat[self.__mean_ren_mask, :])[:, self.__mean_ren_mask]

        w, p = np.linalg.eig(CMat_ren)
        if not np.all(w > 0):
            print 'using pseudo-inverse covariance matrix for "ren" group in region {tag}'.format(tag=tag)
            self.data_ren_covs[tag] = np.linalg.pinv(CMat_ren)
        else:
            self.data_ren_covs[tag] = np.linalg.inv(CMat_ren)

        self.data_convs[tag] = ConvMat


    def __import_theory(self, tables, counterterms, db, my_config):

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

        data = (model, growth_params, loop_params, XY_params, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum)

        with sqlite3.connect(db) as conn:

            for tag in tables:

                ell0 = self.__import_theory_P_ell(conn, tag, data, 0)
                ell2 = self.__import_theory_P_ell(conn, tag, data, 2)
                ell4 = self.__import_theory_P_ell(conn, tag, data, 4)

                self.theory_payload[tag] = np.array([ell0, ell2, ell4])

            for tag in counterterms:

                self.theory_payload[tag] = self.__import_theory_counterterm(conn, tag, data)

            # need f to compute mu^6 counterterm
            self.__import_f(conn, data)

        # finally, construct stochastic counterterms
        ks = self.WiggleZ_conv_ks
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

        self.theory_payload['d1'] = np.array([d1_P0, d1_P2, d1_P4])
        self.theory_payload['d2'] = np.array([d2_P0, d2_P2, d2_P4])
        self.theory_payload['d3'] = np.array([d3_P0, d3_P2, d3_P4])


    def __import_theory_P_ell(self, conn, tag, data, ell):

        model_id, growth_params_id, loop_id, XY_params_id, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum = data

        cursor = conn.cursor()

        table_name = '{tag}_P{ell}'.format(tag=tag, ell=ell)

        cursor.execute(
            "SELECT k_config.k AS k, sample.P1loopSPT_resum AS Pell FROM (SELECT * FROM " + table_name + " WHERE mid=:model AND growth_params=:growth AND loop_params=:loop AND XY_params=:XY AND zid=:zid AND init_Pk_id=:init_Pk AND final_Pk_id=:final_Pk AND IR_cutoff_id=:IR_cutoff AND UV_cutoff_id=:UV_cutoff AND IR_resum_id=:IR_resum) AS sample INNER JOIN k_config ON sample.kid = k_config.id ORDER BY k;",
            {'model': model_id, 'growth': growth_params_id, 'loop': loop_id, 'XY': XY_params_id, 'zid': zid,
             'init_Pk': init_Pk, 'final_Pk': final_Pk, 'IR_cutoff': IR_cutoff, 'UV_cutoff': UV_cutoff,
             'IR_resum': IR_resum})

        ks = []
        Pells = []

        for row in cursor:

            k, Pell = row
            ks.append(k)
            Pells.append(Pell)

        return np.asarray(Pells)


    def __import_theory_counterterm(self, conn, tag, data):

        cursor = conn.cursor()

        model_id, growth_params_id, loop_id, XY_params_id, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum = data

        table_name = 'counterterms_{tag}'.format(tag=tag)

        cursor.execute(
            "SELECT k_config.k AS k, sample.P0_k2_resum AS P0, sample.P2_k2_resum AS P2, sample.P4_k2_resum AS P4 FROM (SELECT * FROM " + table_name + " WHERE mid=:model AND growth_params=:growth AND XY_params=:XY AND zid=:zid AND init_Pk_id=:init_Pk AND final_Pk_id=:final_Pk AND IR_cutoff_id=:IR_cutoff AND UV_cutoff_id=:UV_cutoff AND IR_resum_id=:IR_resum) AS sample INNER JOIN k_config ON sample.kid = k_config.id ORDER BY k;",
            {'model': model_id, 'growth': growth_params_id, 'XY': XY_params_id, 'zid': zid,
             'init_Pk': init_Pk, 'final_Pk': final_Pk, 'IR_cutoff': IR_cutoff, 'UV_cutoff': UV_cutoff,
             'IR_resum': IR_resum})

        ks = []
        P0s = []
        P2s = []
        P4s = []

        for row in cursor:

            k, P0, P2, P4 = row
            ks.append(k)
            P0s.append(P0)
            P2s.append(P2)
            P4s.append(P4)

        return np.array([ np.asarray(P0s), np.asarray(P2s), np.asarray(P4s) ])


    def __import_f(self, conn, data):

        cursor = conn.cursor()

        model_id, growth_params_id, loop_id, XY_params_id, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum = data

        cursor.execute("SELECT f_linear FROM f_factors WHERE zid=:zid;", {'zid': zid})

        flist = cursor.fetchone()

        if flist is None:
            raise LookupError

        self.__f = flist[0]


    def make_coeffs(self, params):

        # extract parameters
        b1_1 = params['b1_1']
        b1_2 = params['b1_2']
        b1_3 = params['b1_3']
        b2_2 = params['b2_2']
        b2_3 = params['b2_3']
        b3 = params['b3']
        bG2_2 = params['bG2_2']
        bG2_3 = params['bG2_3']
        bdG2 = params['bdG2']
        bGamma3 = params['bGamma3']

        # prepare coefficient dictionary
        coeffs = {'nobias': 1.0,
                  'b1_1': b1_1,
                  'b1_2': b1_2,
                  'b1_3': b1_3,
                  'b2_2': b2_2,
                  'b2_3': b2_3,
                  'b3': b3,
                  'bG2_2': bG2_2,
                  'bG2_3': bG2_3,
                  'bdG2': bdG2,
                  'bGamma3': bGamma3,
                  'b1_1_b1_1': b1_1 * b1_1,
                  'b1_2_b1_2': b1_2 * b1_2,
                  'b1_1_b1_2': b1_1 * b1_2,
                  'b1_1_b1_3': b1_1 * b1_3,
                  'b1_1_b2_2': b1_1 * b2_2,
                  'b1_1_b2_3': b1_1 * b2_3,
                  'b1_2_b2_2': b1_2 * b2_2,
                  'b2_2_b2_2': b2_2 * b2_2,
                  'b1_1_b3': b1_1 * b3,
                  'b1_1_bG2_2': b1_1 * bG2_2,
                  'b1_1_bG2_3': b1_1 * bG2_3,
                  'b1_2_bG2_2': b1_2 * bG2_2,
                  'bG2_2_bG2_2': bG2_2 * bG2_2,
                  'b2_2_bG2_2': b2_2 * bG2_2,
                  'b1_1_bdG2': b1_1 * bdG2,
                  'b1_1_bGamma3': b1_1 * bGamma3}

        return coeffs


    def add_counterterms(self, coeffs, values, blinear):

        fb = self.__f / blinear

        c0 = values['c0']
        c2 = values['c2']
        c4 = values['c4']
        d1 = values['d1']
        d2 = values['d2']
        d3 = values['d3']

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        # update coefficient dictionary with counterterm values
        coeffs.update({'c0': c0 * blinear,
                       'c2': c2 * blinear,
                       'c4': c4 * blinear,
                       'c6': c6_ope * blinear,
                       'd1': d1,
                       'd2': d2,
                       'd3': d3})


    def compute(self, block, params, blinear, likes):

        coeffs = self.make_coeffs(params)

        # compute EFT counterterms
        # note this will update coeffs with the values of counterterms during the optimization process,
        # but these will be overwritten with the bestfit values immediately below
        values = self.__compute_EFT_counterterms(coeffs, blinear)
        self.add_counterterms(coeffs, values, blinear)

        # build theoretical P_ell with these counterterms
        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        # sum likelihood over all regions and store
        lik = sum([self.__compute_likelihood(region, P0, P2, P4, 'all') for region in self.data_regions])
        block[likes, 'GLOBAL_LIKE'] = lik

        # store derived EFT parameters for output with the rest of the chain
        block['counterterms', 'c0'] = values['c0']
        block['counterterms', 'c2'] = values['c2']
        block['counterterms', 'c4'] = values['c4']

        block['counterterms', 'd1'] = values['d1']
        block['counterterms', 'd2'] = values['d2']
        block['counterterms', 'd3'] = values['d3']

        return 0


    def build_theory_P_ell(self, coeffs):

        # construct P_ell, including mixing counterterms and stochastic counterterms
        zip = [ coeffs[key] * data for key, data in self.theory_payload.iteritems() ]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a

        return P[0], P[1], P[2]


    def __compute_likelihood(self, region, P0, P2, P4, type='all'):

        if type is 'all':
            mask = self.conv_all_mask
            means = self.data_all_means[region]
            cov = self.data_all_covs[region]
        else:
            mask = self.__conv_ren_mask
            means = self.data_ren_means[region]
            cov = self.data_ren_covs[region]

        conv = self.data_convs[region]

        # first, convolve theory vector with WiggleZ convolution matrix for this region
        P = np.concatenate( (P0, P2, P4) )
        Pconv = np.dot(conv, P)

        # strip out components that can be compared with measurement
        Ptheory = Pconv[mask]

        # compute chi^2
        Delta = Ptheory - means

        return -np.dot(np.dot(Delta, cov), Delta) / 2.0


    def __compute_EFT_counterterms(self, coeffs, blinear):

        fb = self.__f / blinear

        # optimize fit for counterterms over the renormalization region
        # coeffs is updated with the values of the counterterms
        initial_counterterms = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        res = optimize.minimize(self.__compute_renormalization_fit, initial_counterterms, method='Powell',
                                args=(coeffs, blinear, fb), options={'maxiter': 50000, 'maxfev': 50000})

        if not res.success:
            raise RuntimeError(res.message)

        return {'c0': res.x[0],
                'c2': res.x[1],
                'c4': res.x[2],
                'd1': res.x[3],
                'd2': res.x[4],
                'd3': res.x[5]}


    def __compute_renormalization_fit(self, x, coeffs, blinear, fb):

        c0, c2, c4, d1, d2, d3 = x

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        coeffs.update({'c0': c0 * blinear,
                       'c2': c2 * blinear,
                       'c4': c4 * blinear,
                       'c6': c6_ope * blinear,
                       'd1': d1,
                       'd2': d2,
                       'd3': d3})

        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        lik = sum([self.__compute_likelihood(region, P0, P2, P4, 'ren') for region in self.data_regions])

        return -lik


    def cleanup(self):

        return 0
