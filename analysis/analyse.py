from astropy.io import ascii
import csv
import os
import datetime
from getdist import mcsamples as mcs
from getdist import plots as gdp


def ensure_folder(path):

    # generate folder hierarchy if it does not already exist
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise

    return path


class analyse_core(object):

    def __init__(self, t, p):

        self.__params = p
        self.__tools = t


    def get_bias_params(self):

        return self.__params


    def get_model_params(self):

        return self.__tools.theory.parameters


    def get_tools(self):

        return self.__tools


    def compute_chisquare(self, ps, mp, glb):

        # compute chi-square for the parameters in ps

        param_dict = mp(ps)
        coeffs = self.__tools.make_coeff_dict(param_dict)
        P0, P2, P4 = self.__tools.theory.build_theory_P_ell(coeffs, ps, glb(ps))

        chisq = -2.0 * self.__tools.compute_likelihood(P0, P2, P4, 'fit')

        return chisq


class analyse_cosmosis(analyse_core):

    def __init__(self, t, p, root_path, cosmosis_file, out_file, ensemble_file, mp, glb, mixing_plot=None, stochastic_plot=None):

        super(analyse_cosmosis, self).__init__(t, p)

        # construct hierarchy of plot folders
        mixing_folder = ensure_folder(os.path.join(root_path, 'plots', 'mixing'))
        stochastic_folder = ensure_folder(os.path.join(root_path, 'plots', 'stochastic'))
        getdist_folder = ensure_folder(os.path.join(root_path, 'GetDist_output'))

        # generate paths for GetDist-format output files
        getdist_root_file = os.path.join(getdist_folder, out_file)
        getdist_param_file = os.path.join(getdist_folder, out_file + '.paramnames')
        getdist_chain_file = os.path.join(getdist_folder, out_file + '.txt')

        ensemble_root_file = os.path.join(getdist_folder, ensemble_file)
        ensemble_chain_file = os.path.join(getdist_folder, ensemble_file + '.txt')
        ensemble_lock_file = ensemble_chain_file + '.lock'

        # generate paths for output plots
        triangle_mixing_file = os.path.join(mixing_folder, out_file + '.png')
        triangle_stochastic_file = os.path.join(stochastic_folder, out_file + '.png')


        # CONVERT COSMOSIS FILES TO GETDIST FORMAT

        # write GetDist-compatible .params file, if it does not already exist
        self._write_getdist_params_file(getdist_param_file)

        # generate GetDist-compatible chain file from CosmoSIS output, if it does not already exist
        self.__convert_getdist_chain_file(cosmosis_file, getdist_chain_file)


        # IMPORT CONVERTED CHAIN FILES USING GETDIST

        # import chains files using GetDist and cache its MCSamples object internally
        self.__samples, self.bestfit, self.bestfit_chisquare = self.__bestfit_params(getdist_root_file, mp, glb)

        # use lockfile semaphore to avoid race condition where we read an incompletely-converted chain
        if os.path.exists(ensemble_chain_file) and not os.path.exists(ensemble_lock_file):
            self.__ensemble_samples, self.ensemble_bestfit, self.ensemble_bestfit_chisquare = \
                self.__bestfit_params(ensemble_root_file, mp, glb)
        else:
            self.__ensemble_samples = None
            self.ensemble_bestfit = None
            self.ensemble_bestfit_chisquare = None


        # generate triangle plots of mixing and stochastic parameters if requested
        self.__generate_plots(mixing_plot, stochastic_plot, triangle_mixing_file, triangle_stochastic_file)


    def _write_getdist_params_file(self, file):

        # generate GetDist .paramnames file if it does not already exist
        if not os.path.exists(file):

            # cache list (really OrderedDict) of derived model parameters
            model_params = self.get_model_params()
            bias_params = self.get_bias_params()

            print ':: generating GetDist .paramnames file "{f}"'.format(f=file)

            with open(file, 'w') as g:

                writer = csv.DictWriter(g, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

                # write out parameter definitions for bias model parameters
                for par in bias_params:
                    writer.writerow({'name': par, 'LaTeX': bias_params[par]})

                # write out parameter definitions for derived model parameters (eg. c0 or sigmav)
                for par in model_params:
                    writer.writerow({'name': par + '*', 'LaTeX': model_params[par]})

        else:

            print ':: GetDist .paramnames file "{f}" already exists: leaving intact'.format(f=file)


    def __convert_getdist_chain_file(self, cosmosis_file, getdist_file):

        if not os.path.exists(getdist_file):

            # cache list (really OrderedDict) of derived model parameters
            model_params = self.get_model_params()
            bias_params = self.get_bias_params()

            params = list(bias_params.keys()) + list(model_params.keys())

            print ':: converting CosmoSIS-format chain file "{s}" to GetDist-format chain file "{o}"'.format(s=cosmosis_file, o=getdist_file)

            # note p.keys() must return a list that is ordered in the correct way
            input_columns = params + ['like']
            output_columns = ['weight', 'like'] + params

            table = ascii.read(cosmosis_file, Reader=ascii.NoHeader, names=input_columns)

            # use 'lock' file as a semaphore to avoid race conditions with reading/writing the converted file
            lockfile = getdist_file + '.lock'

            with open(lockfile, 'w') as g:

                g.write(datetime.datetime.now().strftime("%d-%B-%Y %H:%M:%S"))

            with open(getdist_file, 'w') as g:

                writer = csv.DictWriter(g, output_columns, delimiter='\t', restval='MISSING')

                for row in table:

                    row_dict = {}

                    # populate with entries for bias model parameters
                    for par in bias_params:
                        row_dict.update({par: row[par]})

                    # populate with entries for derived model parameters (eg. c0 or sigmav)
                    for par in model_params:
                        row_dict.update({par: row[par]})

                    minus_log_lik = -row['like']
                    row_dict.update({'weight': 1, 'like': minus_log_lik})

                    writer.writerow(row_dict)

            os.remove(lockfile)

        else:

            print ':: GetDist-format chain file "{o}" already exists: leaving intact; no conversion of "{s}"'.format(s=cosmosis_file, o=getdist_file)


    def __bestfit_params(self, getdist_root, mp, glb):

        # cache list (really OrderedDict) of derived model parameters
        model_params = self.get_model_params()
        bias_params = self.get_bias_params()

        params = list(bias_params.keys()) + list(model_params.keys())

        # drop front 40% of rows as burn-in; possibly too conservative
        analysis_settings = {'ignore_rows': 0.4}
        samples = mcs.loadMCSamples(getdist_root, settings=analysis_settings)

        x = samples.getLikeStats()

        r = {p: x.parWithName(p).bestfit_sample for p in params}
        chisquare = self.compute_chisquare(r, mp, glb)

        return samples, r, chisquare


    def __generate_plots(self, mixing_plot, stochastic_plot, mixing_file, stochastic_file):

        if mixing_plot is not None:

            print ':: generating triangle plot for mixing counterterms'

            g = gdp.getSubplotPlotter()
            g.triangle_plot(self.__samples, mixing_plot, shaded=True)
            g.export(mixing_file)

        if stochastic_plot is not None:

            print ':: generating triangle plot for stochastic counterterms'

            h = gdp.getSubplotPlotter()
            h.triangle_plot(self.__samples, stochastic_plot, shaded=True)
            h.export(stochastic_file)


    def get_fit_point(self):

        return self.bestfit


    def get_ensemble_fit_point(self):

        return self.ensemble_bestfit


class analyse_maxlike(analyse_core):

    def __init__(self, t, p, root_path, maxlike_file, ensemble_file, mp, glb):

        super(analyse_maxlike, self).__init__(t, p)

        # construct hierarchy of plot folders
        ensure_folder(os.path.join(root_path, 'plots'))

        # READ OUTPUT FILE GENERATED BY MAXLIKE

        # best-fit for this realization
        self.bestfit, self.best_chisquare = self.__bestfit_params(maxlike_file, mp, glb)

        # best-fit for ensemble average, if present
        if os.path.exists(ensemble_file):
            self.ensemble_bestfit, self.ensemble_best_chisquare = self.__bestfit_params(ensemble_file, mp, glb)
        else:
            self.ensemble_bestfit = None
            self.ensemble_best_chisquare = None


    def __bestfit_params(self, file, mp, glb):

        # cache OrderedDict of model parameters
        model_params = self.get_model_params()
        bias_params = self.get_bias_params()

        params = list(bias_params.keys()) + list(model_params.keys())

        # note p.keys() must return a list that is ordered in the correct way
        input_columns = params + ['like']

        # read output file generated by maxlike sampler.
        # this will consist of a single line giving the best-fit values of all parameters, including derived ones
        table = ascii.read(file, Reader=ascii.NoHeader, names=input_columns)

        row = table[0]

        r = {p: row[p] for p in params}
        chisquare = self.compute_chisquare(r, mp, glb)

        return r, chisquare


    def get_fit_point(self):

        return self.bestfit


    def get_ensemble_fit_point(self):

        return self.ensemble_bestfit


def write_summary(real_list, root_path, out_file, mp, glb):
    """Generates a GetDist-format summary file for all realizations, including the ensemble average if present.

    :param real_list: list-like object containing analysis objects for each realization
    :param root_path: root of output path
    :param out_file: name of output file (without extension)
    :param mp: make_parameters function
    :param glb: get_linear_bias function
    :return: None
    """

    # filenames for GetDist chain-like output
    getdist_param_file = os.path.join(root_path, 'plots', out_file + '.paramnames')
    getdist_chain_file = os.path.join(root_path, 'plots', out_file + '.txt')

    # parameter lists from all realizations should be the same
    bias_params, model_params = __get_parameter_lists(real_list)
    params = list(bias_params.keys()) + list(model_params.keys())

    with open(getdist_param_file, 'w') as f:

        print ':: generating GetDist-format chain summary .paramnames "{f}"'.format(f=getdist_param_file)

        writer = csv.DictWriter(f, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

        for p in bias_params:
            writer.writerow({'name': p, 'LaTeX': bias_params[p]})

        for p in model_params:
            writer.writerow({'name': p + '*', 'LaTeX': model_params[p]})

        writer.writerow({'name': r'bestfit_003*', 'LaTeX': r'\chi^2_{0.03}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_005*', 'LaTeX': r'\chi^2_{0.05}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_007*', 'LaTeX': r'\chi^2_{0.07}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_009*', 'LaTeX': r'\chi^2_{0.09}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_011*', 'LaTeX': r'\chi^2_{0.11}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_013*', 'LaTeX': r'\chi^2_{0.13}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_015*', 'LaTeX': r'\chi^2_{0.15}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_017*', 'LaTeX': r'\chi^2_{0.17}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_029*', 'LaTeX': r'\chi^2_{0.29}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_021*', 'LaTeX': r'\chi^2_{0.21}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_023*', 'LaTeX': r'\chi^2_{0.23}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_025*', 'LaTeX': r'\chi^2_{0.25}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_027*', 'LaTeX': r'\chi^2_{0.27}(\mathrm{best})'})
        writer.writerow({'name': r'bestfit_029*', 'LaTeX': r'\chi^2_{0.29}(\mathrm{best})'})

        writer.writerow({'name': r'ensemble_003*', 'LaTeX': r'\chi^2_{0.03}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_005*', 'LaTeX': r'\chi^2_{0.05}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_007*', 'LaTeX': r'\chi^2_{0.07}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_009*', 'LaTeX': r'\chi^2_{0.09}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_011*', 'LaTeX': r'\chi^2_{0.11}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_013*', 'LaTeX': r'\chi^2_{0.13}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_015*', 'LaTeX': r'\chi^2_{0.15}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_017*', 'LaTeX': r'\chi^2_{0.17}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_029*', 'LaTeX': r'\chi^2_{0.29}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_021*', 'LaTeX': r'\chi^2_{0.21}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_023*', 'LaTeX': r'\chi^2_{0.23}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_025*', 'LaTeX': r'\chi^2_{0.25}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_027*', 'LaTeX': r'\chi^2_{0.27}(\mathrm{ensemble})'})
        writer.writerow({'name': r'ensemble_029*', 'LaTeX': r'\chi^2_{0.29}(\mathrm{ensemble})'})

    with open(getdist_chain_file, 'w') as f:

        print ':: generating GetDist-format chain summary file "{f}"'.format(f=getdist_chain_file)

        chisq_points = ['0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25', '0.27', '0.29']

        best_chisq_points = ['best_' + label for label in chisq_points]
        ensemble_chisq_points = ['ensemble_' + label for label in chisq_points]

        columns = ['weight', 'like'] + params + best_chisq_points + ensemble_chisq_points

        writer = csv.DictWriter(f, columns, delimiter='\t', restval='MISSING')

        for real in real_list:

            # get analysis object for realization, and extract its toolkit
            rlz = real_list[real]
            tools = rlz.get_tools()

            # initial entries in the row dictionary can be the best-fit point
            row = rlz.get_fit_point()

            # turn these parameters into a parameter dictionary and then convert to a coefficient dictionary
            param_dict = mp(row)
            coeffs = tools.make_coeff_dict(param_dict)

            # build P_ell at best-fit-point
            P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, row, glb(row))

            # build cumulative chi-square using best-fit point
            best_deviations = tools.compute_chisq_variation(P0, P2, P4, 'best')
            row.update(best_deviations)

            # get ensemble-average best-fit point, if one exists
            ensemble_best_fit = rlz.get_ensemble_fit_point()

            if ensemble_best_fit is not None:

                # build parameter and coefficient dictionaries
                param_dict = mp(ensemble_best_fit)
                coeffs = tools.make_coeff_dict(param_dict)

                # build P_ell at ensemble-average best-fit point
                P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, ensemble_best_fit, glb(ensemble_best_fit))

                # build cumulative chi-square at ensemble-average best-fit point
                ensemble_deviations = tools.compute_chisq_variation(P0, P2, P4, 'ensemble')
                row.update(ensemble_deviations)

            row.update({'weight': 1, 'like': 1})
            writer.writerow(row)


def __get_parameter_lists(analysis_list):

    bias_params = None
    model_params = None

    for r in analysis_list:

        bp = analysis_list[r].get_bias_params()
        mp = analysis_list[r].get_model_params()

        if bias_params is not None:
            if bias_params != bp:
                raise RuntimeError
        else:
            bias_params = bp

        if model_params is not None:
            if model_params != mp:
                raise RuntimeError
        else:
            model_params = mp

    return bias_params, model_params


def write_Pell(real_list, root_path, out_file, mp, glb):
    """For each realization, write a detailed 'Pell' file containing the fitted power spectrum and other comparisons.

    :param real_list: list-like object containing analysis objects for each realization
    :param root_path: root of output path
    :param out_file: name of output file
    :param mp: make_parameters function object
    :param glb: get_linear_bias function object
    :return: None
    """

    for real in real_list:

        rlz = real_list[real]

        p = os.path.join(root_path, 'plots', out_file + '_' + real + '_Pell.csv')

        bestfit = rlz.get_fit_point()
        tools = rlz.get_tools()

        params = mp(bestfit)
        coeffs = tools.make_coeff_dict(params)
        P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, bestfit, glb(bestfit))

        tools.make_summary_plot(p, P0, P2, P4)
