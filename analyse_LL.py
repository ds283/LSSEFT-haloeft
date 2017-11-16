import numpy as np
from astropy.io import ascii
from scipy import optimize
import sqlite3
import csv
from os import path, makedirs
from getdist import mcsamples as mcs
from getdist import plots as gdp


h1_means = "data/means/pkpole_wizcola_1hr_z0pt2_0pt6.dat"
h3_means = "data/means/pkpole_wizcola_3hr_z0pt2_0pt6.dat"
h9_means = "data/means/pkpole_wizcola_9hr_z0pt2_0pt6.dat"
h11_means = "data/means/pkpole_wizcola_11hr_z0pt2_0pt6.dat"
h15_means = "data/means/pkpole_wizcola_15hr_z0pt2_0pt6.dat"
h22_means = "data/means/pkpole_wizcola_22hr_z0pt2_0pt6.dat"
h1_matrix = "data/covariances/covar_1hr_z0pt2_0pt6.dat"
h3_matrix = "data/covariances/covar_3hr_z0pt2_0pt6.dat"
h9_matrix = "data/covariances/covar_9hr_z0pt2_0pt6.dat"
h11_matrix = "data/covariances/covar_11hr_z0pt2_0pt6.dat"
h15_matrix = "data/covariances/covar_15hr_z0pt2_0pt6.dat"
h22_matrix = "data/covariances/covar_22hr_z0pt2_0pt6.dat"
realization = 1
theory_db = "theory/WizCOLA_CAMB_halo_LL@z=0_kWiggleZ.sqlite"
model = 0
growth_params = 0
loop_params = 0
XY_params = 0
zid = 0
init_Pk = 0
final_Pk = 1
IR_cutoff = 0
UV_cutoff = 0
IR_resum = 0
renormalize_kmin = 0.14
renormalize_kmax = 0.26


class haloeft:

    def __init__(self):

        self.__nbin = 15
        self.__nbinc = 25

        self.__WiggleZ_mean_ks = np.linspace(0.01, 0.29, self.__nbin)
        self.__WiggleZ_conv_ks = np.linspace(0.01, 0.49, self.__nbinc)


        # READ AND CACHE THEORY

        self.__theory_payload = {}

        self.__import_theory(['nobias', 'b1', 'b2', 'bs2', 'b3nl',
                              'b1_b1', 'b2_b2', 'bs2_bs2',
                              'b1_b2', 'b1_bs2', 'b1_b3nl', 'b2_bs2'],
                             ['c0', 'c2', 'c4', 'c6'], theory_db)

        mean_all_mask = self.__WiggleZ_mean_ks > 0.01
        conv_all_mask = np.all([self.__WiggleZ_conv_ks > 0.01, self.__WiggleZ_conv_ks < 0.30], axis=0)
        self.__mean_all_mask = np.concatenate((mean_all_mask, mean_all_mask, mean_all_mask))
        self.__conv_all_mask = np.concatenate((conv_all_mask, conv_all_mask, conv_all_mask))

        mean_ren_mask = np.all([self.__WiggleZ_mean_ks > renormalize_kmin, self.__WiggleZ_mean_ks < renormalize_kmax], axis=0)
        conv_ren_mask = np.all([self.__WiggleZ_conv_ks > renormalize_kmin, self.__WiggleZ_conv_ks < renormalize_kmax], axis=0)
        self.__mean_ren_mask = np.concatenate((mean_ren_mask, mean_ren_mask, mean_ren_mask))
        self.__conv_ren_mask = np.concatenate((conv_ren_mask, conv_ren_mask, conv_ren_mask))


        # READ AND CACHE WIGGLEZ DATA PRODUCTS

        self.__data_raw_means = {}
        self.__data_all_means = {}
        self.__data_ren_means = {}
        self.__data_all_covs = {}
        self.__data_ren_covs = {}
        self.__data_convs = {}
        self.__data_regions = ['1h', '3h', '8h', '11h', '15h', '22h']

        self.__import_WiggleZ_data("1h", h1_means, h1_matrix)
        self.__import_WiggleZ_data("3h", h3_means, h3_matrix)
        self.__import_WiggleZ_data("8h", h9_means, h9_matrix)
        self.__import_WiggleZ_data("11h", h11_means, h11_matrix)
        self.__import_WiggleZ_data("15h", h15_means, h15_matrix)
        self.__import_WiggleZ_data("22h", h22_means, h22_matrix)


    def __import_WiggleZ_data(self, tag, means_file, matrix_file):

        # number of bins in P_ell measurements and covariance matrix
        nbin = self.__nbin

        # number of bins in convolution  matrix
        nbinc = self.__nbinc

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

        self.__data_raw_means[tag] = np.concatenate( (this_P0, this_P2, this_P4) )
        self.__data_all_means[tag] = np.concatenate( (this_P0, this_P2, this_P4) )[self.__mean_all_mask]
        self.__data_ren_means[tag] = np.concatenate( (this_P0, this_P2, this_P4) )[self.__mean_ren_mask]

        CMat = np.empty((3*nbin, 3*nbin))
        ConvMat = np.empty((3*nbinc, 3*nbinc))

        for row in cov_table:

            CMat[row['i']-1, row['j']-1] = row['value']

        for row in conv_table:

            ConvMat[row['i']-1, row['j']-1] = row['value']

        CMat_all = (CMat[self.__mean_all_mask, :])[:, self.__mean_all_mask]

        w, p = np.linalg.eig(CMat_all)
        if not np.all(w > 0):
            print 'using pseudo-inverse covariance matrix for "all" group in region {tag}'.format(tag=tag)
            self.__data_all_covs[tag] = np.linalg.pinv(CMat_all)
        else:
            self.__data_all_covs[tag] = np.linalg.inv(CMat_all)

        CMat_ren = (CMat[self.__mean_ren_mask, :])[:, self.__mean_ren_mask]

        w, p = np.linalg.eig(CMat_ren)
        if not np.all(w > 0):
            print 'using pseudo-inverse covariance matrix for "ren" group in region {tag}'.format(tag=tag)
            self.__data_ren_covs[tag] = np.linalg.pinv(CMat_ren)
        else:
            self.__data_ren_covs[tag] = np.linalg.inv(CMat_ren)

        self.__data_convs[tag] = ConvMat


    def __import_theory(self, tables, counterterms, db):

        data = (model, growth_params, loop_params, XY_params, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum)

        with sqlite3.connect(db) as conn:

            for tag in tables:

                ell0 = self.__import_theory_P_ell(conn, tag, data, 0)
                ell2 = self.__import_theory_P_ell(conn, tag, data, 2)
                ell4 = self.__import_theory_P_ell(conn, tag, data, 4)

                self.__theory_payload[tag] = (ell0, ell2, ell4)

            for tag in counterterms:

                self.__theory_payload[tag] = self.__import_theory_counterterm(conn, tag, data)

            # need f to compute mu^6 counterterm
            self.__import_f(conn, data)

        # finally, construct stochastic counterterms
        ks = self.__WiggleZ_conv_ks
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

        self.__theory_payload['d1'] = (d1_P0, d1_P2, d1_P4)
        self.__theory_payload['d2'] = (d2_P0, d2_P2, d2_P4)
        self.__theory_payload['d3'] = (d3_P0, d3_P2, d3_P4)


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

        return np.asarray(P0s), np.asarray(P2s), np.asarray(P4s)


    def __import_f(self, conn, data):

        cursor = conn.cursor()

        model_id, growth_params_id, loop_id, XY_params_id, zid, init_Pk, final_Pk, IR_cutoff, UV_cutoff, IR_resum = data

        cursor.execute("SELECT f_linear FROM f_factors WHERE zid=:zid;", {'zid': zid})

        flist = cursor.fetchone()

        if flist is None:
            raise LookupError

        self.__f = flist[0]
        print 'f = {f}'.format(f=self.__f)


    def compute_chisq(self, b1, b2, c0, c2, c4, d1, d2, d3):

        fb = self.__f / b1

        bs2 = -4.0*(b1-1.0)/7.0
        b3nl = 32.0*(b1-1.0)/315.0

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        coeffs = {'nobias': 1.0, 'b1': b1, 'b2': b2, 'bs2': bs2, 'b3nl': b3nl,
                  'b1_b1': b1 * b1, 'b2_b2': b2 * b2, 'bs2_bs2': bs2 * bs2,
                  'b1_b2': b1 * b2, 'b1_bs2': b1 * bs2, 'b1_b3nl': b1 * b3nl, 'b2_bs2': b2 * bs2,
                  'c0': c0 * b1, 'c2': c2 * b1, 'c4': c4 * b1, 'c6': c6_ope * b1,
                  'd1': d1, 'd2': d2, 'd3': d3}

        P0, P2, P4 = self.__build_theory_P_ell(coeffs)

        P = np.concatenate( (P0, P2, P4) )

        rval = { }
        labels = ['0.03', '0.05', '0.07', '0.09',
                  '0.11', '0.13', '0.15', '0.17', '0.19',
                  '0.21', '0.23', '0.25', '0.27', '0.29']

        i_mask = np.array([ False for i in xrange(np.sum(self.__mean_all_mask)/3) ])

        for i in xrange(len(i_mask)):

            i_mask[i] = True
            i_mask_full = np.concatenate( (i_mask, i_mask, i_mask) )

            chisq = 0.0

            for region in self.__data_regions:

                mask = self.__conv_all_mask
                means = self.__data_all_means[region]
                cov = self.__data_all_covs[region]
                conv = self.__data_convs[region]

                Pconv = np.dot(conv, P)

                Ptheory = Pconv[mask]

                Ptheory_cut = Ptheory[i_mask_full]
                means_cut = means[i_mask_full]
                cov_cut = cov[i_mask_full,:][:,i_mask_full]

                Delta_cut = Ptheory_cut - means_cut

                chisq += np.dot(np.dot(Delta_cut, cov_cut), Delta_cut)

            rval[labels[i]] = chisq

        return rval


    def make_plot(self, p, b1, b2, c0, c2, c4, d1, d2, d3):

        fb = self.__f / b1

        bs2 = -4.0*(b1-1.0)/7.0
        b3nl = 32.0*(b1-1.0)/315.0

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        coeffs = {'nobias': 1.0, 'b1': b1, 'b2': b2, 'bs2': bs2, 'b3nl': b3nl,
                  'b1_b1': b1 * b1, 'b2_b2': b2 * b2, 'bs2_bs2': bs2 * bs2,
                  'b1_b2': b1 * b2, 'b1_bs2': b1 * bs2, 'b1_b3nl': b1 * b3nl, 'b2_bs2': b2 * bs2,
                  'c0': c0 * b1, 'c2': c2 * b1, 'c4': c4 * b1, 'c6': c6_ope * b1,
                  'd1': d1, 'd2': d2, 'd3': d3}

        P0, P2, P4 = self.__build_theory_P_ell(coeffs)
        P = np.concatenate( (P0, P2, P4) )

        with open(p, 'w') as f:

            labels = ['k', 'P0', 'P2', 'P4']

            for r in self.__data_regions:
                labels = labels + [r + '_P0_theory', r + '_P2_theory', r + '_P4_theory',
                                   r + '_P0_WizCOLA', r + '_P2_WizCOLA', r + '_P4_WizCOLA',
                                   r + '_P0_Delta', r + '_P2_Delta', r + '_P4_Delta']

            writer = csv.DictWriter(f, labels)
            writer.writeheader()

            for i in xrange(len(self.__WiggleZ_mean_ks)):

                row = {'k': self.__WiggleZ_mean_ks[i], 'P0': P0[i], 'P2': P2[i], 'P4': P4[i]}

                for r in self.__data_regions:

                    means = self.__data_raw_means[r]
                    conv = self.__data_convs[r]

                    Pconv = np.dot(conv, P)

                    row.update({r+'_P0_theory': Pconv[i], r+'_P0_WizCOLA': means[i]})
                    row.update({r+'_P2_theory': Pconv[1*self.__nbinc+i], r+'_P2_WizCOLA': means[1*self.__nbin+i]})
                    row.update({r+'_P4_theory': Pconv[2*self.__nbinc+i], r+'_P4_WizCOLA': means[2*self.__nbin+i]})

                    row.update({r+'_P0_Delta': row[r+'_P0_theory'] - row[r+'_P0_WizCOLA']})
                    row.update({r+'_P2_Delta': row[r+'_P2_theory'] - row[r+'_P2_WizCOLA']})
                    row.update({r+'_P4_Delta': row[r+'_P4_theory'] - row[r+'_P4_WizCOLA']})

                writer.writerow(row)


    def __build_theory_P_ell(self, coeffs):

        # construct P_ell, including mixing counterterms and stochastic counterterms
        P0 = np.zeros( (self.__nbinc,) )
        P2 = np.zeros( (self.__nbinc,) )
        P4 = np.zeros( (self.__nbinc,) )

        for table, data in self.__theory_payload.iteritems():

            P0 = P0 + coeffs[table] * data[0]
            P2 = P2 + coeffs[table] * data[1]
            P4 = P4 + coeffs[table] * data[2]

        return P0, P2, P4


    def __compute_lik(self, region, P0, P2, P4, type='all'):

        if type is 'all':
            mask = self.__conv_all_mask
            means = self.__data_all_means[region]
            cov = self.__data_all_covs[region]
        else:
            mask = self.__conv_ren_mask
            means = self.__data_ren_means[region]
            cov = self.__data_ren_covs[region]

        conv = self.__data_convs[region]

        # first, convolve theory vector with WiggleZ convolution matrix for this region
        P = np.concatenate( (P0, P2, P4) )
        Pconv = np.dot(conv, P)

        # strip out components that can be compared with measurement
        Ptheory = Pconv[mask]

        # compute chi^2
        Delta = Ptheory - means

        return -np.dot(np.dot(Delta, cov), Delta) / 2.0


    def __compute_EFT_counterterms(self, b1, b2, bs2, b3nl, fb):

        # optimize fit for counterterms over the renormalization region

        # prepare coefficient dictionary without counterterms; these will be filled in by the
        # __compute_EFT_likelihood function and overwritten on each step
        coeffs = {'nobias': 1.0, 'b1': b1, 'b2': b2, 'bs2': bs2, 'b3nl': b3nl,
                  'b1_b1': b1 * b1, 'b2_b2': b2 * b2, 'bs2_bs2': bs2 * bs2,
                  'b1_b2': b1 * b2, 'b1_bs2': b1 * bs2, 'b1_b3nl': b1 * b3nl, 'b2_bs2': b2 * bs2}

        res = optimize.minimize(self.__compute_EFT_likelihood, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], method='Powell',
                                args=(coeffs, b1, fb), options={'maxiter': 14000, 'maxfev': 14000})

        if not res.success:
            raise RuntimeError(res.message)

        return res.x


    def __compute_EFT_likelihood(self, x, coeffs, b1, fb):

        c0, c2, c4, d1, d2, d3 = x

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        EFT_coeffs = {'c0': c0 * b1, 'c2': c2 * b1, 'c4': c4 * b1, 'c6': c6_ope * b1, 'd1': d1, 'd2': d2, 'd3': d3}
        coeffs.update(EFT_coeffs)

        P0, P2, P4 = self.__build_theory_P_ell(coeffs)

        lik = sum([self.__compute_lik(region, P0, P2, P4, 'ren') for region in self.__data_regions])

        return -lik


class analyse:

    def __init__(self, emcee_root, root):

        emcee_file = path.join('output', emcee_root)

        table = ascii.read(emcee_file, Reader=ascii.NoHeader,
                           names=['b1', 'b2', 'c0', 'c2', 'c4', 'd1', 'd2', 'd3', 'like'])

        mixing_folder = path.join('plots', 'mixing')
        stochastic_folder = path.join('plots', 'stochastic')

        if not path.exists(mixing_folder):
            makedirs(mixing_folder)

        if not path.exists(stochastic_folder):
            makedirs(stochastic_folder)

        param_names_file = path.join('plots', root + '.paramnames')
        getdist_root = path.join('plots', root)
        getdist_file = path.join('plots', root + '.txt')
        triangle_mixing_file = path.join(mixing_folder, root + '.png')
        triangle_stochastic_file = path.join(stochastic_folder, root + '.png')

        # generate GetDist .paramnames file
        if not path.exists(param_names_file):

            with open(param_names_file, 'w') as g:

                writer = csv.DictWriter(g, ['name', 'LaTeX'], delimiter='\t')

                writer.writerow({'name': r'b1', 'LaTeX': r'b_1'})
                writer.writerow({'name': r'b2', 'LaTeX': r'b_2'})
                writer.writerow({'name': r'c0*', 'LaTeX': r'c_0'})
                writer.writerow({'name': r'c2*', 'LaTeX': r'c_2'})
                writer.writerow({'name': r'c4*', 'LaTeX': r'c_4'})
                writer.writerow({'name': r'd1*', 'LaTeX': r'd_1'})
                writer.writerow({'name': r'd2*', 'LaTeX': r'd_2'})
                writer.writerow({'name': r'd3*', 'LaTeX': r'd_3'})

        if not path.exists(getdist_file):

            with open(getdist_file, 'w') as g:

                writer = csv.DictWriter(g, ['weight', 'like', 'b1', 'b2', 'c0', 'c2', 'c4', 'd1', 'd2', 'd3'], delimiter='\t')

                for row in table:

                    writer.writerow({'weight': 1, 'like': -row['like'], 'b1': row['b1'], 'b2': row['b2'],
                                     'c0': row['c0'], 'c2': row['c2'], 'c4': row['c4'],
                                     'd1': row['d1'], 'd2': row['d2'], 'd3': row['d3']})


        # import chains files using GetDist and cache its MCSamples object internally
        analysis_settings = {'ignore_rows': 0.4}
        self.__samples = mcs.loadMCSamples(getdist_root, settings=analysis_settings)

        g = gdp.getSubplotPlotter()
        g.triangle_plot(self.__samples, ['b1', 'b2', 'c0', 'c2', 'c4'], shaded=True)
        g.export(triangle_mixing_file)

        h = gdp.getSubplotPlotter()
        h.triangle_plot(self.__samples, ['b1', 'b2', 'd1', 'd2', 'd3'], shaded=True)
        h.export(triangle_stochastic_file)


    def get_fit_point(self):

        x = self.__samples.getLikeStats()

        return (x.parWithName('b1').bestfit_sample,
                x.parWithName('b2').bestfit_sample,
                x.parWithName('c0').bestfit_sample,
                x.parWithName('c2').bestfit_sample,
                x.parWithName('c4').bestfit_sample,
                x.parWithName('d1').bestfit_sample,
                x.parWithName('d2').bestfit_sample,
                x.parWithName('d3').bestfit_sample)



def write_summary(list, halo_model, root):

    # filenames for GetDist chain-like output
    param_file = path.join('plots', root + '.paramnames')
    summary_file = path.join('plots', root + '.txt')

    with open(param_file, 'w') as f:

        writer = csv.DictWriter(f, ['name', 'LaTeX'], delimiter='\t')

        writer.writerow({'name': r'b1', 'LaTeX': r'b_1'})
        writer.writerow({'name': r'b2', 'LaTeX': r'b_2'})
        writer.writerow({'name': r'c0*', 'LaTeX': r'c_0'})
        writer.writerow({'name': r'c2*', 'LaTeX': r'c_2'})
        writer.writerow({'name': r'c4*', 'LaTeX': r'c_4'})
        writer.writerow({'name': r'd1*', 'LaTeX': r'd_1'})
        writer.writerow({'name': r'd2*', 'LaTeX': r'd_2'})
        writer.writerow({'name': r'd3*', 'LaTeX': r'd_3'})
        writer.writerow({'name': r'chisq003*', 'LaTeX': r'\chi^2_{0.03}'})
        writer.writerow({'name': r'chisq005*', 'LaTeX': r'\chi^2_{0.05}'})
        writer.writerow({'name': r'chisq007*', 'LaTeX': r'\chi^2_{0.07}'})
        writer.writerow({'name': r'chisq009*', 'LaTeX': r'\chi^2_{0.09}'})
        writer.writerow({'name': r'chisq011*', 'LaTeX': r'\chi^2_{0.11}'})
        writer.writerow({'name': r'chisq013*', 'LaTeX': r'\chi^2_{0.13}'})
        writer.writerow({'name': r'chisq015*', 'LaTeX': r'\chi^2_{0.15}'})
        writer.writerow({'name': r'chisq017*', 'LaTeX': r'\chi^2_{0.17}'})
        writer.writerow({'name': r'chisq029*', 'LaTeX': r'\chi^2_{0.29}'})
        writer.writerow({'name': r'chisq021*', 'LaTeX': r'\chi^2_{0.21}'})
        writer.writerow({'name': r'chisq023*', 'LaTeX': r'\chi^2_{0.23}'})
        writer.writerow({'name': r'chisq025*', 'LaTeX': r'\chi^2_{0.25}'})
        writer.writerow({'name': r'chisq027*', 'LaTeX': r'\chi^2_{0.27}'})
        writer.writerow({'name': r'chisq029*', 'LaTeX': r'\chi^2_{0.29}'})

    with open(summary_file, 'w') as f:

        writer = csv.DictWriter(f, ['weight', 'like',
                                    'b1', 'b2', 'c0', 'c2', 'c4', 'd1', 'd2', 'd3',
                                    '0.03', '0.05', '0.07', '0.09',
                                    '0.11', '0.13', '0.15', '0.17', '0.19',
                                    '0.21', '0.23', '0.25', '0.27', '0.29'], delimiter='\t')

        for real in list:

            b1, b2, c0, c2, c4, d1, d2, d3 = list[real].get_fit_point()

            row = {'weight': 1, 'like': 1,
                   'b1': b1, 'b2': b2, 'c0': c0, 'c2': c2, 'c4': c4, 'd1': d1, 'd2': d2, 'd3': d3}

            deviations = halo_model.compute_chisq(b1, b2, c0, c2, c4, d1, d2, d3)
            row.update(deviations)

            writer.writerow(row)


def write_Pell(list, halo_model, root):

    for real in list:

        p = path.join('plots', root + '_' + real + '_Pell.csv')

        b1, b2, c0, c2, c4, d1, d2, d3 = list[real].get_fit_point()

        halo_model.make_plot(p, b1, b2, c0, c2, c4, d1, d2, d3)



halo_model = haloeft()

r01 = analyse('output_LL_r01.txt', 'LL_r01')
r02 = analyse('output_LL_r02.txt', 'LL_r02')
r03 = analyse('output_LL_r03.txt', 'LL_r03')
r04 = analyse('output_LL_r04.txt', 'LL_r04')
r05 = analyse('output_LL_r05.txt', 'LL_r05')
r06 = analyse('output_LL_r06.txt', 'LL_r06')
r07 = analyse('output_LL_r07.txt', 'LL_r07')
r08 = analyse('output_LL_r08.txt', 'LL_r08')
r09 = analyse('output_LL_r09.txt', 'LL_r09')
r10 = analyse('output_LL_r10.txt', 'LL_r10')

list = {'r01': r01, 'r02': r02, 'r03': r03, 'r04': r04, 'r05': r05, 'r06': r06, 'r07': r07, 'r08': r08, 'r09': r09, 'r10': r10}

write_summary(list, halo_model, 'LL_ensemble')
write_Pell(list, halo_model, 'LL_ensemble')
