import numpy as np
from astropy.io import ascii
import csv
import os
from getdist import mcsamples as mcs
from getdist import plots as gdp

import haloeft as heft

import WizCOLA
import EFT

PconvCeiling=1E6


class EFT_tools(heft.cosmosis_pipeline):

    def __init__(self, realization, model_name, fit_kmin=0.02, fit_kmax=0.30, ren_kmin=0.02, ren_kmax=0.30):

        self.realization = realization
        self.model_name = model_name

        # use a dictionary to mimic the CosmoSIS datablock API
        config = {}
        
        root = os.path.join("..", "..")

        config["HaloEFT", "h1_means"] = os.path.join(root, "data/means/pkpole_wizcola_1hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h3_means"] = os.path.join(root, "data/means/pkpole_wizcola_3hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h9_means"] = os.path.join(root, "data/means/pkpole_wizcola_9hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h11_means"] = os.path.join(root, "data/means/pkpole_wizcola_11hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h15_means"] = os.path.join(root, "data/means/pkpole_wizcola_15hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h22_means"] = os.path.join(root, "data/means/pkpole_wizcola_22hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h1_matrix"] = os.path.join(root, "data/covariances/covar_1hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h3_matrix"] = os.path.join(root, "data/covariances/covar_3hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h9_matrix"] = os.path.join(root, "data/covariances/covar_9hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h11_matrix"] = os.path.join(root, "data/covariances/covar_11hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h15_matrix"] = os.path.join(root, "data/covariances/covar_15hr_z0pt2_0pt6.dat")
        config["HaloEFT", "h22_matrix"] = os.path.join(root, "data/covariances/covar_22hr_z0pt2_0pt6.dat")
        config["HaloEFT", "realization"] = realization
        config["HaloEFT", "theory_db"] = os.path.join(root, "theory/WizCOLA_CAMB_halo_full@z=0_kWiggleZ.sqlite")
        config["HaloEFT", "model"] = 0
        config["HaloEFT", "growth_params"] = 0
        config["HaloEFT", "loop_params"] = 0
        config["HaloEFT", "XY_params"] = 0
        config["HaloEFT", "zid"] = 0
        config["HaloEFT", "init_Pk"] = 0
        config["HaloEFT", "final_Pk"] = 1
        config["HaloEFT", "IR_cutoff"] = 0
        config["HaloEFT", "UV_cutoff"] = 0
        config["HaloEFT", "IR_resum"] = 0
        config["HaloEFT", "fit_kmin"] = fit_kmin
        config["HaloEFT", "fit_kmax"] = fit_kmax
        config["HaloEFT", "renormalize_kmin"] = ren_kmin
        config["HaloEFT", "renormalize_kmax"] = ren_kmax

        # set up container of k-samples for WizCOLA
        ks = WizCOLA.ksamples(config)

        # build data container
        data = WizCOLA.products(config, ks)

        # build theory container
        theory = EFT.theory(config, ks)

        # pass configuration to base class
        super(EFT_tools, self).__init__('EFT_Analyse', data, theory)


    def compute_chisq_variation(self, P0, P2, P4):

        P = np.concatenate( (P0, P2, P4) )

        # set up empty dictionary to hold return values
        rval = {}

        # set up mask: we will gradually add True values to this as we step through the different
        # k-sample points that contribute to the final chi-square value
        # notice the mask is for just one P_ell individually, not the concatenated group (P0, P2, P4)
        i_mask = np.array([False for i in xrange(len(self.data.k_sample.WiggleZ_mean_ks))])

        convolved_Pk = {}
        # cache convolved theory power spectra
        for r in self.data.regions:

            Pconv = np.dot(self.data.convs[r], P)

            if np.any(abs(Pconv) > PconvCeiling):
                print '!! region {tag}: Pconv very large: largest value = {l}, P.max = {P}, conv.max = {c}'.format(
                    tag=r, l=abs(Pconv).max(), P=abs(P).max(), c=abs(self.data.convs[r]).max())
                print Pconv
                raise RuntimeError

            convolved_Pk[r] = Pconv[self.data.k_sample.conv_to_means_mask]

        # step through each available sample point
        for i in xrange(len(i_mask)):

            # is this sample point included in the fit mask?
            i_mask[i] = self.data.k_sample.mean_fit_mask[i]

            # if so, compute the chi-square up to this k-sample point; otherwise, ignore it
            if(i_mask[i]):

                # extend mask for the concatenated group (P0, P2, P4)
                i_mask_full = np.concatenate( (i_mask, i_mask, i_mask) )

                # zero accumulator for chi-square
                chisq = 0.0

                # loop over all data regions:
                for r in self.data.regions:

                    # extract data products for this region
                    means = self.data.fit_means[r]
                    inv_cov = self.data.fit_inv_covs[r]

                    # cut raw region used in fit (ie. from 0.01 to 0.29) out of full convolved vector
                    Ptheory = convolved_Pk[r]

                    # cut down to just the k-sample points we are using here
                    Ptheory_cut = Ptheory[i_mask_full]

                    # cut a mask for the data means and data covariance matrix out of the mask for the raw region
                    i_mask_means = i_mask_full[self.data.k_sample.mean_fit_mask]
                    means_cut = means[i_mask_means]
                    inv_cov_cut = inv_cov[i_mask_means,:][:,i_mask_means]

                    # take difference between theory and data, and process into a chi-square
                    Delta_cut = Ptheory_cut - means_cut
                    chisq += np.dot(np.dot(Delta_cut, inv_cov_cut), Delta_cut)

                rval[self.data.k_sample.labels[i]] = chisq

        return rval


    def make_summary_plot(self, plot_file, P0, P2, P4):

        P = np.concatenate( (P0, P2, P4) )

        # compute total chi-square value
        total_chisq = -2.0 * self.compute_likelihood(P0, P2, P4, 'fit')

        with open(plot_file, 'w') as f:

            print ':: generating theory + data fit comparison file "{f}"'.format(f=plot_file)

            # add labels for global columns
            labels = ['ireal', 'chisq', 'model',
                      'k', 'P0', 'P2', 'P4',
                      'P0_theory_avg', 'P2_theory_avg', 'P4_theory_avg',
                      'P0_WizCOLA_avg', 'P2_WizCOLA_avg', 'P4_WizCOLA_avg',
                      'P0_WizCOLA_err_avg', 'P2_WizCOLA_err_avg', 'P4_WizCOLA_err_avg',
                      'P0_ensemble_avg', 'P2_ensemble_avg', 'P4_ensemble_avg',
                      'P0_ensemble_err_avg', 'P2_ensemble_err_avg', 'P4_ensemble_err_avg',
                      'P0_Delta_avg', 'P2_Delta_avg', 'P4_Delta_avg']

            # supplement with labels for per-region columns
            for r in self.data.regions:

                labels += [r + '_P0_theory', r + '_P2_theory', r + '_P4_theory',
                           r + '_P0_WizCOLA', r + '_P2_WizCOLA', r + '_P4_WizCOLA',
                           r + '_P0_WizCOLA_err', r + '_P2_WizCOLA_err', r + '_P4_WizCOLA_err',
                           r + '_P0_Delta', r + '_P2_Delta', r + '_P4_Delta']

            # open a CSV writer for these columns and emit a header
            writer = csv.DictWriter(f, labels, restval='MISSING')
            writer.writeheader()

            convolved_Pk = {}
            # cache convolved theory power spectra
            for r in self.data.regions:

                Pconv = np.dot(self.data.convs[r], P)

                if np.any(abs(Pconv) > PconvCeiling):
                    print '!! region {tag}: Pconv very large: largest value = {l}, P.max = {P}, conv.max = {c}'.format(
                        tag=r, l=abs(Pconv).max(), P=abs(P).max(), c=abs(self.data.convs[r]).max())
                    print Pconv
                    raise RuntimeError

                convolved_Pk[r] = Pconv

            # loop over each k sample point to be included in the analysis,
            # and then process its contribution for each region
            for i in xrange(len(self.data.k_sample.WiggleZ_mean_ks)):

                row = {'k': self.data.k_sample.WiggleZ_mean_ks[i], 'P0': P0[i], 'P2': P2[i], 'P4': P4[i]}

                if i == 0:
                    row.update({'ireal': self.realization, 'chisq': total_chisq, 'model': self.model_name})
                else:
                    row.update({'ireal': -99, 'chisq': -99, 'model': '-'})

                # initialize accumulators used to sum over all regions
                P0_theory_tot = 0
                P2_theory_tot = 0
                P4_theory_tot = 0

                P0_WizCOLA_tot = 0
                P2_WizCOLA_tot = 0
                P4_WizCOLA_tot = 0

                P0_ensemble_tot = 0
                P2_ensemble_tot = 0
                P4_ensemble_tot = 0

                P0_WizCOLA_var_tot = 0
                P2_WizCOLA_var_tot = 0
                P4_WizCOLA_var_tot = 0

                P0_ensemble_var_tot = 0
                P2_ensemble_var_tot = 0
                P4_ensemble_var_tot = 0

                P0_Delta_tot = 0
                P2_Delta_tot = 0
                P4_Delta_tot = 0

                # loop over regions
                for r in self.data.regions:

                    means = self.data.raw_means[r]
                    variances = self.data.raw_variance[r]

                    ensemble_means = self.data.ensemble_means[r]
                    ensemble_variances = self.data.ensemble_variance[r]

                    # read previously cached convolved power spectra
                    Pconv = convolved_Pk[r]

                    # pick out P0, P2, P4 values for this k-sample point
                    P0_theory = Pconv[0 * self.data.k_sample.nbinc + i]
                    P2_theory = Pconv[1 * self.data.k_sample.nbinc + i]
                    P4_theory = Pconv[2 * self.data.k_sample.nbinc + i]

                    P0_data = means[0 * self.data.k_sample.nbin + i]
                    P2_data = means[1 * self.data.k_sample.nbin + i]
                    P4_data = means[2 * self.data.k_sample.nbin + i]

                    P0_ensemble = ensemble_means[0 * self.data.k_sample.nbin + i]
                    P2_ensemble = ensemble_means[1 * self.data.k_sample.nbin + i]
                    P4_ensemble = ensemble_means[2 * self.data.k_sample.nbin + i]

                    P0_data_var = variances[0 * self.data.k_sample.nbin + i]
                    P2_data_var = variances[1 * self.data.k_sample.nbin + i]
                    P4_data_var = variances[2 * self.data.k_sample.nbin + i]

                    P0_ensemble_var = ensemble_variances[0 * self.data.k_sample.nbin + i]
                    P2_ensemble_var = ensemble_variances[1 * self.data.k_sample.nbin + i]
                    P4_ensemble_var = ensemble_variances[2 * self.data.k_sample.nbin + i]

                    row.update({r + '_P0_theory': P0_theory, r + '_P0_WizCOLA': P0_data, r + '_P0_WizCOLA_err': np.sqrt(P0_data_var)})
                    row.update({r + '_P2_theory': P2_theory, r + '_P2_WizCOLA': P2_data, r + '_P2_WizCOLA_err': np.sqrt(P2_data_var)})
                    row.update({r + '_P4_theory': P4_theory, r + '_P4_WizCOLA': P4_data, r + '_P4_WizCOLA_err': np.sqrt(P4_data_var)})

                    P0_delta = P0_data - P0_theory
                    P2_delta = P2_data - P2_theory
                    P4_delta = P4_data - P4_theory

                    row.update({r + '_P0_Delta': P0_delta})
                    row.update({r + '_P2_Delta': P2_delta})
                    row.update({r + '_P4_Delta': P4_delta})

                    P0_theory_tot += P0_theory
                    P2_theory_tot += P2_theory
                    P4_theory_tot += P4_theory

                    P0_WizCOLA_tot += P0_data
                    P2_WizCOLA_tot += P2_data
                    P4_WizCOLA_tot += P4_data

                    P0_ensemble_tot += P0_ensemble
                    P2_ensemble_tot += P2_ensemble
                    P4_ensemble_tot += P4_ensemble

                    P0_WizCOLA_var_tot += P0_data_var
                    P2_WizCOLA_var_tot += P2_data_var
                    P4_WizCOLA_var_tot += P4_data_var

                    P0_ensemble_var_tot += P0_ensemble_var
                    P2_ensemble_var_tot += P2_ensemble_var
                    P4_ensemble_var_tot += P4_ensemble_var

                    P0_Delta_tot += P0_delta
                    P2_Delta_tot += P2_delta
                    P4_Delta_tot += P4_delta

                regions = len(self.data.regions)

                row.update({'P0_theory_avg': P0_theory_tot/regions,
                            'P2_theory_avg': P2_theory_tot/regions,
                            'P4_theory_avg': P4_theory_tot/regions,
                            'P0_WizCOLA_avg': P0_WizCOLA_tot/regions,
                            'P2_WizCOLA_avg': P2_WizCOLA_tot/regions,
                            'P4_WizCOLA_avg': P4_WizCOLA_tot/regions,
                            'P0_WizCOLA_err_avg': np.sqrt(P0_WizCOLA_var_tot)/regions,
                            'P2_WizCOLA_err_avg': np.sqrt(P2_WizCOLA_var_tot)/regions,
                            'P4_WizCOLA_err_avg': np.sqrt(P4_WizCOLA_var_tot)/regions,
                            'P0_ensemble_avg': P0_ensemble_tot/regions,
                            'P2_ensemble_avg': P2_ensemble_tot/regions,
                            'P4_ensemble_avg': P4_ensemble_tot/regions,
                            'P0_ensemble_err_avg': np.sqrt(P0_ensemble_var_tot)/regions,
                            'P2_ensemble_err_avg': np.sqrt(P2_ensemble_var_tot)/regions,
                            'P4_ensemble_err_avg': np.sqrt(P4_ensemble_var_tot)/regions,
                            'P0_Delta_avg': P0_Delta_tot/regions,
                            'P2_Delta_avg': P2_Delta_tot/regions,
                            'P4_Delta_avg': P4_Delta_tot/regions})

                writer.writerow(row)


class analyse_core(object):

    def __init__(self, r, mn, p, mp, glb):

        self.realization = r
        self.params = p

        self.make_params = mp
        self.get_linear_bias = glb

        self.tools = EFT_tools(r, mn)


    def get_params(self):

        return self.params


    def get_tools(self):

        return self.tools


    def compare_chisquare(self, ps):

        # compute chi-square for the parameters in ps

        tools = self.tools

        param_dict = self.make_params(ps)
        coeffs = tools.make_coeff_dict(param_dict)
        P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, ps, self.get_linear_bias(ps))

        chisq = -2.0 * tools.compute_likelihood(P0, P2, P4, 'fit')

        return chisq


class analyse_CosmoSIS(analyse_core):

    def __init__(self, r, mn, p, CosmoSIS_file, out_file, make_params, get_linear_bias, mixing_plot=None, stochastic_plot=None):

        super(analyse_CosmoSIS, self).__init__(r, mn, p, make_params, get_linear_bias)

        # construct hierarchy of plot folders
        mixing_folder = os.path.join('plots', 'mixing')
        stochastic_folder = os.path.join('plots', 'stochastic')

        # generate folder hierarchy if it does not already exist
        if not os.path.exists(mixing_folder):
            try:
                os.makedirs(mixing_folder)
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise

        if not os.path.exists(stochastic_folder):
            try:
                os.makedirs(stochastic_folder)
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise

        # generate paths for GetDist-format output files
        getdist_root = os.path.join('plots', out_file)
        getdist_param_file = os.path.join('plots', out_file + '.paramnames')
        getdist_chain_file = os.path.join('plots', out_file + '.txt')

        # generate paths for output plots
        triangle_mixing_file = os.path.join(mixing_folder, out_file + '.png')
        triangle_stochastic_file = os.path.join(stochastic_folder, out_file + '.png')


        # CONVERT COSMOSIS FILES TO GETDIST FORMAT

        # generate GetDist .paramnames file if it does not already exist
        if not os.path.exists(getdist_param_file):

            print ':: generating GetDist .paramnames file "{f}"'.format(f=getdist_param_file)

            with open(getdist_param_file, 'w') as g:

                writer = csv.DictWriter(g, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

                for par in p:
                    writer.writerow({'name': par, 'LaTeX': p[par]})

                writer.writerow({'name': r'c0*', 'LaTeX': r'c_0'})
                writer.writerow({'name': r'c2*', 'LaTeX': r'c_2'})
                writer.writerow({'name': r'c4*', 'LaTeX': r'c_4'})
                writer.writerow({'name': r'd1*', 'LaTeX': r'd_1'})
                writer.writerow({'name': r'd2*', 'LaTeX': r'd_2'})
                writer.writerow({'name': r'd3*', 'LaTeX': r'd_3'})

        else:

            print ':: GetDist .paramnames file "{f}" already exists: leaving intact'.format(f=getdist_param_file)


        # generate GetDist-compatible chain file from CosmoSIS output, if it does not already exist
        CosmoSIS_path = os.path.join('output', CosmoSIS_file)

        if not os.path.exists(getdist_chain_file):

            print ':: converting CosmoSIS-format chain file "{s}" to GetDist-format chain file "{o}"'.format(s=CosmoSIS_path, o=getdist_chain_file)

            # note p.keys() must return a list that is ordered in the correct way
            input_columns = list(p.keys()) + ['c0', 'c2', 'c4', 'd1', 'd2', 'd3', 'like']
            output_columns = ['weight', 'like'] + list(p.keys()) + ['c0', 'c2', 'c4', 'd1', 'd2', 'd3']

            table = ascii.read(CosmoSIS_path, Reader=ascii.NoHeader, names=input_columns)

            with open(getdist_chain_file, 'w') as g:

                writer = csv.DictWriter(g, output_columns, delimiter='\t', restval='MISSING')

                for row in table:

                    row_dict = {'c0': row['c0'], 'c2': row['c2'], 'c4': row['c4'],
                                'd1': row['d1'], 'd2': row['d2'], 'd3': row['d3']}

                    for par in p:
                        row_dict.update({par: row[par]})

                    minus_log_lik = -row['like']
                    row_dict.update({'weight': 1, 'like': minus_log_lik})

                    writer.writerow(row_dict)

        else:

            print ':: GetDist-format chain file "{o}" already exists: leaving intact; no conversion of "{s}"'.format(s=CosmoSIS_path, o=getdist_chain_file)


        # IMPORT CONVERTED CHAIN FILES USING GETDIST

        # import chains files using GetDist and cache its MCSamples object internally
        analysis_settings = {'ignore_rows': 0.4}
        self.__samples = mcs.loadMCSamples(getdist_root, settings=analysis_settings)

        if mixing_plot is not None:

            print ':: generating triangle plot for mixing counterterms'

            g = gdp.getSubplotPlotter()
            g.triangle_plot(self.__samples, mixing_plot, shaded=True)
            g.export(triangle_mixing_file)

        if stochastic_plot is not None:

            print ':: generating triangle plot for stochastic counterterms'

            h = gdp.getSubplotPlotter()
            h.triangle_plot(self.__samples, stochastic_plot, shaded=True)
            h.export(triangle_stochastic_file)

        x = self.__samples.getLikeStats()

        r = {p: x.parWithName(p).bestfit_sample for p in self.params}

        EFT_r = {'c0': x.parWithName('c0').bestfit_sample,
                 'c2': x.parWithName('c2').bestfit_sample,
                 'c4': x.parWithName('c4').bestfit_sample,
                 'd1': x.parWithName('d1').bestfit_sample,
                 'd2': x.parWithName('d2').bestfit_sample,
                 'd3': x.parWithName('d3').bestfit_sample}

        r.update(EFT_r)

        self.bestfit = r
        self.bestfit_chisquare = self.compare_chisquare(r)


    def get_fit_point(self):

        return self.bestfit


class analyse_maxlike(analyse_core):

    def __init__(self, r, mn, p, maxlike_file, make_params, get_linear_bias):

        super(analyse_maxlike, self).__init__(r, mn, p, make_params, get_linear_bias)

        # construct hierarchy of plot folders
        plot_folder = os.path.join('plots')

        if not os.path.exists(plot_folder):
            try:
                os.makedirs(plot_folder)
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise


        # READ OUTPUT FILE GENERATED BY MAXLIKE

        # this will consist of a single line giving the best-fit values of all parameters, including derived ones
        maxlike_path = os.path.join('output', maxlike_file)

        # note p.keys() must return a list that is ordered in the correct way
        input_columns = list(p.keys()) + ['c0', 'c2', 'c4', 'd1', 'd2', 'd3', 'like']

        table = ascii.read(maxlike_path, Reader=ascii.NoHeader, names=input_columns)

        row = table[0]

        r = {p: row[p] for p in list(p.keys() + ['c0', 'c2', 'c4', 'd1', 'd2', 'd3'])}

        self.bestfit = r
        self.best_chisquare = self.compare_chisquare(r)


    def get_fit_point(self):

        return self.bestfit



def write_summary(realizations, out_file):

    # filenames for GetDist chain-like output
    getdist_param_file = os.path.join('plots', out_file + '.paramnames')
    getdist_chain_file = os.path.join('plots', out_file + '.txt')

    # parameter list from all realizations should be the same
    params = None
    for r in realizations:
        p = realizations[r].get_params()

        if params is not None:
            if params != p:
                raise RuntimeError
        else:
            params = p

    with open(getdist_param_file, 'w') as f:

        print ':: generating GetDist-format chain summary .paramnames "{f}"'.format(f=getdist_param_file)

        writer = csv.DictWriter(f, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

        for p in params:
            writer.writerow({'name': p, 'LaTeX': params[p]})

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

    with open(getdist_chain_file, 'w') as f:

        print ':: generating GetDist-format chain summary file "{f}"'.format(f=getdist_chain_file)

        columns = ['weight', 'like'] + list(params.keys()) + \
                  ['c0', 'c2', 'c4', 'd1', 'd2', 'd3'] + \
                  ['0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25', '0.27', '0.29']

        writer = csv.DictWriter(f, columns, delimiter='\t', restval='MISSING')

        for real in realizations:

            rlz = realizations[real]

            row = rlz.get_fit_point()
            tools = rlz.get_tools()

            param_dict = rlz.make_params(row)
            coeffs = tools.make_coeff_dict(param_dict)
            P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, row, rlz.get_linear_bias(row))

            deviations = tools.compute_chisq_variation(P0, P2, P4)

            row.update(deviations)
            row.update({'weight': 1, 'like': 1})

            writer.writerow(row)


def write_Pell(list, out_file):

    for real in list:

        rlz = list[real]

        p = os.path.join('plots', out_file + '_' + real + '_Pell.csv')

        bestfit = rlz.get_fit_point()
        tools = rlz.get_tools()

        params = rlz.make_params(bestfit)
        coeffs = tools.make_coeff_dict(params)
        P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, bestfit, rlz.get_linear_bias(bestfit))

        tools.make_summary_plot(p, P0, P2, P4)
