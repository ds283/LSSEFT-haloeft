import numpy as np

from astropy.io import ascii
from scipy import optimize

import sqlite3

import settings
import utilities as utils


# container class for k-sample points
class WiggleZ_ksamples(object):

    def __init__(self, my_config):

        self.nbin = 15
        self.nbinc = 25

        self.WiggleZ_mean_ks = np.linspace(0.01, 0.29, self.nbin)
        self.WiggleZ_conv_ks = np.linspace(0.01, 0.49, self.nbinc)

        self.labels = ['0.01', '0.03', '0.05', '0.07', '0.09',
                       '0.11', '0.13', '0.15', '0.17', '0.19',
                       '0.21', '0.23', '0.25', '0.27', '0.29']


        # BUILD AND CACHE MASKS TO SELECT DIFFERENT PARTS OF THE DATA PRODUCTS

        fit_kmin = my_config["HaloEFT", "fit_kmin"]
        fit_kmax = my_config["HaloEFT", "fit_kmax"]

        mam = np.all([self.WiggleZ_mean_ks > fit_kmin, self.WiggleZ_mean_ks <= fit_kmax], axis=0)
        cam = np.all([self.WiggleZ_conv_ks > fit_kmin, self.WiggleZ_conv_ks <= fit_kmax], axis=0)
        self.mean_fit_mask = np.concatenate((mam, mam, mam))
        self.conv_fit_mask = np.concatenate((cam, cam, cam))

        ren_kmin = my_config["HaloEFT", "renormalize_kmin"]
        ren_kmax = my_config["HaloEFT", "renormalize_kmax"]

        mrm = np.all([self.WiggleZ_mean_ks > ren_kmin, self.WiggleZ_mean_ks <= ren_kmax], axis=0)
        crm = np.all([self.WiggleZ_conv_ks > ren_kmin, self.WiggleZ_conv_ks <= ren_kmax], axis=0)
        self.mean_ren_mask = np.concatenate((mrm, mrm, mrm))
        self.conv_ren_mask = np.concatenate((crm, crm, crm))

        cmm = self.WiggleZ_conv_ks <= 0.30
        self.conv_to_means_mask = np.concatenate((cmm, cmm, cmm))


# container class for WizCOLA data products
class WizCOLA_products(object):

    def __init__(self, my_config, k_sample):

        self.k_sample = k_sample


        # READ AND CACHE WIGGLEZ DATA PRODUCTS

        h_means = {}
        h_cov = {}

        h_means.update({'1h': my_config["HaloEFT", "h1_means"]})
        h_cov.update({'1h': my_config["HaloEFT", "h1_matrix"]})

        h_means.update({'3h': my_config["HaloEFT", "h3_means"]})
        h_cov.update({'3h': my_config["HaloEFT", "h3_matrix"]})

        h_means.update({'9h': my_config["HaloEFT", "h9_means"]})
        h_cov.update({'9h': my_config["HaloEFT", "h9_matrix"]})

        h_means.update({'11h': my_config["HaloEFT", "h11_means"]})
        h_cov.update({'11h': my_config["HaloEFT", "h11_matrix"]})

        h_means.update({'15h': my_config["HaloEFT", "h15_means"]})
        h_cov.update({'15h': my_config["HaloEFT", "h15_matrix"]})

        h_means.update({'22h': my_config["HaloEFT", "h22_means"]})
        h_cov.update({'22h': my_config["HaloEFT", "h22_matrix"]})

        # object-local storage for data products, held as dictionaries indexed by region label

        # 1. Power spectrum mean values: (a) raw, (b) cut down to region needed for likelihood fit,
        # (c) cut down to region needed for renormalization
        self.raw_means = {}
        self.fit_means = {}
        self.ren_means = {}

        # 2. Inverse covariance matrices: (a) cut down to region needed for likelihood fit,
        # (b) cut down to region needed for renormalization
        self.fit_inv_covs = {}
        self.ren_inv_covs = {}

        # 4. Measurement errors on the power spectrum mean value
        self.raw_variance = {}

        # 5. Ensemble average measurements, per region
        self.ensemble_means = {}
        self.ensemble_variance = {}

        # 6. Convolution matrix taking theory outputs to each WiggleZ region
        # (eg. accounts for geometry and selection function)
        self.convs = {}

        # 7. Region labels
        self.regions = ['1h', '3h', '9h', '11h', '15h', '22h']

        # perform import
        for r in self.regions:

            self.__import_WiggleZ_data(r, h_means[r], h_cov[r], my_config)


    def __import_WiggleZ_data(self, tag, means_file, matrix_file, my_config):

        # extract desired realization from datablock
        realization = my_config["HaloEFT", "realization"]

        # number of bins in P_ell measurements and covariance matrix
        nbin = self.k_sample.nbin

        # number of bins in convolution  matrix
        nbinc = self.k_sample.nbinc

        # read in the WiggleZ data for this region as astropy.table.Table instances
        mean_table = ascii.read(means_file, Reader=ascii.NoHeader,
                                names=['ireal', 'k', 'P0', 'P0err', 'P2', 'P2err', 'P4', 'P4err'])

        cov_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=0, data_end=(3*nbin)**2,
                               names=['i', 'j', 'value'])

        conv_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=(3*nbin)**2, data_end=(3*nbin)**2+(3*nbinc)**2,
                                names=['i', 'j', 'value'])

        # want to extract the power spectrum for a given realization, so group mean_table by realization
        mean_by_ireal = mean_table.group_by('ireal')

        # get list of all realizations
        all_realizations = mean_by_ireal.groups.keys['ireal']

        # use mask to extract Pk means for just this realization
        mean_this_realization = mean_by_ireal.groups[realization]

        # sort into ascending order of k
        mean_this_realization.sort('k')

        # extract P0, P2, P4 columns as numpy arrays
        this_P0 = np.asarray(mean_this_realization['P0'])
        this_P2 = np.asarray(mean_this_realization['P2'])
        this_P4 = np.asarray(mean_this_realization['P4'])

        # concatenate these arrays to form a single block consisting of P0 values followed by P2, P4 values
        Pblock = np.concatenate((this_P0, this_P2, this_P4))

        # extract relevant components -- 'raw' is everything, but 'fit' and 'ren' use only a subset of values
        self.raw_means[tag] = Pblock
        self.fit_means[tag] = Pblock[self.k_sample.mean_fit_mask]
        self.ren_means[tag] = Pblock[self.k_sample.mean_ren_mask]


        # construct empty matrices to hold covariance and convolution matrices after read-in
        covariance = np.empty((3*nbin, 3*nbin))
        convolution = np.empty((3*nbinc, 3*nbinc))

        # import matrices
        utils.import_matrix(covariance, cov_table, tag, 'covariance')
        utils.import_matrix(convolution, conv_table, tag, 'convolution')


        # strip out diagonal entries of covariance matrix as error estimates
        variances = np.diagonal(covariance)
        self.raw_variance[tag] = variances


        # cut covariance matrix down only to those values used in the likelihood fit
        # (notice we do this *before* inversion)
        covariance_cut_fit = (covariance[self.k_sample.mean_fit_mask, :])[:, self.k_sample.mean_fit_mask]

        # invert and store
        self.fit_inv_covs[tag] = utils.invert_covariance_matrix(covariance_cut_fit, tag, 'all')


        # cut covariance matrix down only to those values used in renormalization
        # (notice we do this *before* inversion)
        covariance_cut_ren = (covariance[self.k_sample.mean_ren_mask, :])[:, self.k_sample.mean_ren_mask]

        # invert and store
        self.ren_inv_covs[tag] = utils.invert_covariance_matrix(covariance_cut_ren, tag, 'ren')

        # store convolution matrix, first checking whether it contains any anomalously-large elements
        # that would be symptomatic of an error
        if np.any(abs(convolution) > settings.ConvCeiling):

            print '!! region {tag}: convolution matrix contains very large element = {v}'.format(tag=tag, v=abs(convolution).max())
            raise RuntimeError

        self.convs[tag] = convolution


        # extract ensemble average power spectrum for this region, which is stored as realization 0
        ensemble_means = mean_by_ireal.groups[0]
        ensemble_means.sort('k')

        ensemble_P0 = np.asarray(ensemble_means['P0'])
        ensemble_P2 = np.asarray(ensemble_means['P2'])
        ensemble_P4 = np.asarray(ensemble_means['P4'])

        ensemble_Pk = np.concatenate((ensemble_P0, ensemble_P2, ensemble_P4))

        self.ensemble_means[tag] = ensemble_Pk
        self.ensemble_variance[tag] = variances / (len(all_realizations) - 1.0)


# container class for EFT theory products
class EFT_products(object):

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


class HaloEFT_core(object):

    def __init__(self, my_config, my_name, data, theory):

        self.__mod_name = my_name

        self.data = data
        self.theory = theory


    def make_coeff_dict(self, params):

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


    def add_counterterms_to_coeff_dict(self, coeffs, values, blinear):

        fb = self.theory.f / blinear

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

        coeffs = self.make_coeff_dict(params)

        # compute EFT counterterms
        # note this will update coeffs with the values of counterterms during the optimization process,
        # but these will be overwritten with the bestfit values immediately below
        values = self.__compute_EFT_counterterms(coeffs, blinear)
        self.add_counterterms_to_coeff_dict(coeffs, values, blinear)

        # build theoretical P_ell with these counterterms
        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        # sum likelihood over all regions and store back into the datablock
        lik = sum([self.compute_likelihood(region, P0, P2, P4, 'fit') for region in self.data.regions])
        block[likes, 'HALOEFT_LIKE'] = lik


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

        # this implementation is a bit cryptic but has been tuned for speed -- this is the rate-limiting
        # step in an MCMC analysis
        zip = [ coeffs[key] * data for key, data in self.theory.payload.iteritems() ]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a              # in-place addition is fastest

        return P[0], P[1], P[2]


    def compute_likelihood(self, region, P0, P2, P4, type='fit'):

        if type is 'fit':
            mask = self.data.k_sample.conv_fit_mask
            means = self.data.fit_means[region]
            inv_cov = self.data.fit_inv_covs[region]
        else:
            mask = self.data.k_sample.conv_ren_mask
            means = self.data.ren_means[region]
            inv_cov = self.data.ren_inv_covs[region]

        conv = self.data.convs[region]

        # first, convolve theory vector with WiggleZ convolution matrix for this region
        P = np.concatenate( (P0, P2, P4) )
        Pconv = np.dot(conv, P)

        # strip out components that can be compared with measurement
        Ptheory = Pconv[mask]

        # compute chi^2
        Delta = Ptheory - means
        return -np.dot(np.dot(Delta, inv_cov), Delta) / 2.0


    def __compute_EFT_counterterms(self, coeffs, blinear):

        fb = self.theory.f / blinear

        # optimize fit for counterterms over the renormalization region
        # coeffs is updated with the values of the counterterms
        initial_counterterms = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        res = optimize.minimize(self.__compute_renormalization_fit, initial_counterterms, method='Powell',
                                args=(coeffs, blinear, fb), options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 50000, 'maxfev': 50000})

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

        lik = sum([self.compute_likelihood(region, P0, P2, P4, 'ren') for region in self.data.regions])

        return -lik


    def cleanup(self):

        return 0
