import numpy as np
from astropy.io import ascii
import csv
from os import path, makedirs
from getdist import mcsamples as mcs
from getdist import plots as gdp
import haloeft as heft


class EFT_tools(heft.HaloEFT_core):

    def __init__(self, realization):

        # use a dictionary to mimic the CosmoSIS datablock API
        config = {}

        config["HaloEFT", "h1_means"] = "data/means/pkpole_wizcola_1hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h3_means"] = "data/means/pkpole_wizcola_3hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h9_means"] = "data/means/pkpole_wizcola_9hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h11_means"] = "data/means/pkpole_wizcola_11hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h15_means"] = "data/means/pkpole_wizcola_15hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h22_means"] = "data/means/pkpole_wizcola_22hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h1_matrix"] = "data/covariances/covar_1hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h3_matrix"] = "data/covariances/covar_3hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h9_matrix"] = "data/covariances/covar_9hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h11_matrix"] = "data/covariances/covar_11hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h15_matrix"] = "data/covariances/covar_15hr_z0pt2_0pt6.dat"
        config["HaloEFT", "h22_matrix"] = "data/covariances/covar_22hr_z0pt2_0pt6.dat"
        config["HaloEFT", "realization"] = realization
        config["HaloEFT", "theory_db"] = "theory/WizCOLA_CAMB_halo_full@z=0_kWiggleZ.sqlite"
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
        config["HaloEFT", "renormalize_kmin"] = 0.02
        config["HaloEFT", "renormalize_kmax"] = 0.30

        # pass configuration to base class
        # this arranges for read-in of the WiggleZ data products and contents of the theory database
        super(EFT_tools, self).__init__(config, 'EFT_Analyse')


    def compute_chisq(self, coeffs):

        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        P = np.concatenate( (P0, P2, P4) )

        rval = {}
        labels = ['0.03', '0.05', '0.07', '0.09',
                  '0.11', '0.13', '0.15', '0.17', '0.19',
                  '0.21', '0.23', '0.25', '0.27', '0.29']

        i_mask = np.array([False for i in xrange(np.sum(self.mean_all_mask) / 3)])

        for i in xrange(len(i_mask)):

            i_mask[i] = True
            i_mask_full = np.concatenate( (i_mask, i_mask, i_mask) )

            chisq = 0.0

            for region in self.data_regions:

                mask = self.conv_all_mask
                means = self.data_all_means[region]
                cov = self.data_all_covs[region]
                conv = self.data_convs[region]

                Pconv = np.dot(conv, P)

                Ptheory = Pconv[mask]

                Ptheory_cut = Ptheory[i_mask_full]
                means_cut = means[i_mask_full]
                cov_cut = cov[i_mask_full,:][:,i_mask_full]

                Delta_cut = Ptheory_cut - means_cut

                chisq += np.dot(np.dot(Delta_cut, cov_cut), Delta_cut)

            rval[labels[i]] = chisq

        return rval


    def make_plot(self, plot_file, coeffs):

        P0, P2, P4 = self.build_theory_P_ell(coeffs)
        P = np.concatenate( (P0, P2, P4) )

        with open(plot_file, 'w') as f:

            print 'generating theory + fit file "{f}"'.format(f=plot_file)

            labels = ['k', 'P0', 'P2', 'P4',
                      'P0_theory_avg', 'P2_theory_avg', 'P4_theory_avg',
                      'P0_WizCOLA_avg', 'P2_WizCOLA_avg', 'P4_WizCOLA_avg',
                      'P0_Delta_avg', 'P2_Delta_avg', 'P4_Delta_avg']

            for r in self.data_regions:
                labels = labels + [r + '_P0_theory', r + '_P2_theory', r + '_P4_theory',
                                   r + '_P0_WizCOLA', r + '_P2_WizCOLA', r + '_P4_WizCOLA',
                                   r + '_P0_Delta', r + '_P2_Delta', r + '_P4_Delta']

            writer = csv.DictWriter(f, labels)
            writer.writeheader()

            for i in xrange(len(self.WiggleZ_mean_ks)):

                row = {'k': self.WiggleZ_mean_ks[i], 'P0': P0[i], 'P2': P2[i], 'P4': P4[i]}

                P0_theory_tot = 0
                P2_theory_tot = 0
                P4_theory_tot = 0

                P0_WizCOLA_tot = 0
                P2_WizCOLA_tot = 0
                P4_WizCOLA_tot = 0

                P0_Delta_tot = 0
                P2_Delta_tot = 0
                P4_Delta_tot = 0

                for r in self.data_regions:

                    means = self.data_raw_means[r]
                    conv = self.data_convs[r]

                    Pconv = np.dot(conv, P)

                    row.update({r+'_P0_theory': Pconv[i], r+'_P0_WizCOLA': means[i]})
                    row.update({r+'_P2_theory': Pconv[1 * self.nbinc + i], r + '_P2_WizCOLA': means[1 * self.nbin + i]})
                    row.update({r+'_P4_theory': Pconv[2 * self.nbinc + i], r + '_P4_WizCOLA': means[2 * self.nbin + i]})

                    row.update({r+'_P0_Delta': row[r+'_P0_theory'] - row[r+'_P0_WizCOLA']})
                    row.update({r+'_P2_Delta': row[r+'_P2_theory'] - row[r+'_P2_WizCOLA']})
                    row.update({r+'_P4_Delta': row[r+'_P4_theory'] - row[r+'_P4_WizCOLA']})

                    P0_theory_tot += row[r+'_P0_theory']
                    P2_theory_tot += row[r+'_P2_theory']
                    P4_theory_tot += row[r+'_P4_theory']

                    P0_WizCOLA_tot += row[r+'_P0_WizCOLA']
                    P2_WizCOLA_tot += row[r+'_P2_WizCOLA']
                    P4_WizCOLA_tot += row[r+'_P4_WizCOLA']

                    P0_Delta_tot += P0_theory_tot - P0_WizCOLA_tot
                    P2_Delta_tot += P2_theory_tot - P2_WizCOLA_tot
                    P4_Delta_tot += P4_theory_tot - P4_WizCOLA_tot

                regions = len(self.data_regions)

                row.update({'P0_theory_avg': P0_theory_tot/regions,
                            'P2_theory_avg': P2_theory_tot/regions,
                            'P4_theory_avg': P4_theory_tot/regions,
                            'P0_WizCOLA_avg': P0_WizCOLA_tot/regions,
                            'P2_WizCOLA_avg': P2_WizCOLA_tot/regions,
                            'P4_WizCOLA_avg': P4_WizCOLA_tot/regions,
                            'P0_Delta_avg': P0_Delta_tot/regions,
                            'P2_Delta_avg': P2_Delta_tot/regions,
                            'P4_Delta_avg': P4_Delta_tot/regions})

                writer.writerow(row)


class analyse(object):

    def __init__(self, realization, emcee_file, out_file, params, mixing_plot=None, stochastic_plot=None):

        self.__realization = realization
        self.__params = params
        self.__tools = EFT_tools(realization)

        emcee_path = path.join('output', emcee_file)

        input_columns = list(params.keys()) + ['c0', 'c2', 'c4', 'd1', 'd2', 'd3', 'like']

        mixing_folder = path.join('plots', 'mixing')
        stochastic_folder = path.join('plots', 'stochastic')

        if not path.exists(mixing_folder):
            makedirs(mixing_folder)

        if not path.exists(stochastic_folder):
            makedirs(stochastic_folder)

        param_names_file = path.join('plots', out_file + '.paramnames')
        getdist_root = path.join('plots', out_file)
        getdist_file = path.join('plots', out_file + '.txt')
        triangle_mixing_file = path.join(mixing_folder, out_file + '.png')
        triangle_stochastic_file = path.join(stochastic_folder, out_file + '.png')


        # generate GetDist .paramnames file if it does not already exist
        if not path.exists(param_names_file):

            print 'generating GetDist .paramnames file "{f}"'.format(f=param_names_file)

            with open(param_names_file, 'w') as g:

                writer = csv.DictWriter(g, ['name', 'LaTeX'], delimiter='\t')

                for p in params:
                    writer.writerow({'name': p, 'LaTeX': params[p]})

                writer.writerow({'name': r'c0*', 'LaTeX': r'c_0'})
                writer.writerow({'name': r'c2*', 'LaTeX': r'c_2'})
                writer.writerow({'name': r'c4*', 'LaTeX': r'c_4'})
                writer.writerow({'name': r'd1*', 'LaTeX': r'd_1'})
                writer.writerow({'name': r'd2*', 'LaTeX': r'd_2'})
                writer.writerow({'name': r'd3*', 'LaTeX': r'd_3'})

        else:

            print 'GetDist .paramnames file "{f}" already exists: leaving intact'.format(f=param_names_file)


        # generate GetDist-compatible chain file from emcee output, if it does not already exist
        if not path.exists(getdist_file):

            print 'converting emcee chain file "{s}" to GetDist-format chain file "{o}"'.format(s=emcee_path, o=getdist_file)

            table = ascii.read(emcee_path, Reader=ascii.NoHeader, names=input_columns)

            with open(getdist_file, 'w') as g:

                writer = csv.DictWriter(g, ['weight', 'like', 'b1', 'b2', 'c0', 'c2', 'c4', 'd1', 'd2', 'd3'], delimiter='\t')

                for row in table:

                    row_dict = {'weight': 1,
                                'like': -row['like'],
                                'c0': row['c0'],
                                'c2': row['c2'],
                                'c4': row['c4'],
                                'd1': row['d1'],
                                'd2': row['d2'],
                                'd3': row['d3']}

                    for p in params:
                        row_dict.update({p: row[p]})

                    writer.writerow(row_dict)

        else:

            print 'GetDist-format chain file "{o}" already exists: leaving intact; no conversion of "{s}"'.format(s=emcee_path, o=getdist_file)


        # import chains files using GetDist and cache its MCSamples object internally
        analysis_settings = {'ignore_rows': 0.4}
        self.__samples = mcs.loadMCSamples(getdist_root, settings=analysis_settings)

        if mixing_plot is not None:

            print 'generating triangle plot for mixing counterterms'

            g = gdp.getSubplotPlotter()
            g.triangle_plot(self.__samples, mixing_plot, shaded=True)
            g.export(triangle_mixing_file)

        if stochastic_plot is not None:

            print 'generating triangle plot for stochastic counterterms'

            h = gdp.getSubplotPlotter()
            h.triangle_plot(self.__samples, stochastic_plot, shaded=True)
            h.export(triangle_stochastic_file)


    def get_params(self):

        return self.__params


    def get_tools(self):

        return self.__tools


    def get_fit_point(self):

        x = self.__samples.getLikeStats()

        r = {p: x.parWithName(p).bestfit_sample for p in self.__params}

        EFT_r = {'c0': x.parWithName('c0').bestfit_sample,
                 'c2': x.parWithName('c2').bestfit_sample,
                 'c4': x.parWithName('c4').bestfit_sample,
                 'd1': x.parWithName('d1').bestfit_sample,
                 'd2': x.parWithName('d2').bestfit_sample,
                 'd3': x.parWithName('d3').bestfit_sample}

        r.update(EFT_r)

        return r


def write_summary(realizations, make_params, get_blinear, out_file):

    # filenames for GetDist chain-like output
    param_file = path.join('plots', out_file + '.paramnames')
    summary_file = path.join('plots', out_file + '.txt')

    # parameter list from all realizations should be the same
    params = None
    for r in realizations:
        p = realizations[r].get_params()

        if params is not None:
            if params != p:
                raise RuntimeError
        else:
            params = p

    with open(param_file, 'w') as f:

        print 'generating GetDist-format chain summary .paramnames "{f}"'.format(f=param_file)

        writer = csv.DictWriter(f, ['name', 'LaTeX'], delimiter='\t')

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

    with open(summary_file, 'w') as f:

        print 'generating GetDist-format chain summary file "{f}"'.format(f=summary_file)

        columns = ['weight', 'like'] + list(params.keys()) + \
                  ['c0', 'c2', 'c4', 'd1', 'd2', 'd3'] + \
                  ['0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25', '0.27', '0.29']

        writer = csv.DictWriter(f, columns, delimiter='\t')

        for real in realizations:

            row = realizations[real].get_fit_point()
            tools = realizations[real].get_tools()

            param_dict = make_params(row)
            coeffs = tools.make_coeffs(param_dict)
            tools.add_counterterms(coeffs, row, get_blinear(row))

            deviations = tools.compute_chisq(coeffs)

            row.update(deviations)
            row.update({'weight': 1, 'like': 1})

            writer.writerow(row)


def write_Pell(list, make_params, get_blinear, out_file):

    for real in list:

        p = path.join('plots', out_file + '_' + real + '_Pell.csv')

        bestfit = list[real].get_fit_point()
        tools = list[real].get_tools()

        params = make_params(bestfit)
        coeffs = tools.make_coeffs(params)
        tools.add_counterterms(coeffs, bestfit, get_blinear(bestfit))

        tools.make_plot(p, coeffs)
