from astropy.io import ascii
import csv
import os
from getdist import mcsamples as mcs
from getdist import plots as gdp


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

    def __init__(self, t, p, root_path, cosmosis_file, out_file, mp, glb, mixing_plot=None, stochastic_plot=None):

        super(analyse_cosmosis, self).__init__(t, p)

        # construct hierarchy of plot folders
        mixing_folder = os.path.join(root_path, 'plots', 'mixing')
        stochastic_folder = os.path.join(root_path, 'plots', 'stochastic')

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
        getdist_root = os.path.join(root_path, 'plots', out_file)
        getdist_param_file = os.path.join(root_path, 'plots', out_file + '.paramnames')
        getdist_chain_file = os.path.join(root_path, 'plots', out_file + '.txt')

        # generate paths for output plots
        triangle_mixing_file = os.path.join(mixing_folder, out_file + '.png')
        triangle_stochastic_file = os.path.join(stochastic_folder, out_file + '.png')


        # CONVERT COSMOSIS FILES TO GETDIST FORMAT

        # cache list (really OrderedDict) of derived model parameters
        model_params = self.get_model_params()

        # generate GetDist .paramnames file if it does not already exist
        if not os.path.exists(getdist_param_file):

            print ':: generating GetDist .paramnames file "{f}"'.format(f=getdist_param_file)

            with open(getdist_param_file, 'w') as g:

                writer = csv.DictWriter(g, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

                # write out parameter definitions for bias model parameters
                for par in p:
                    writer.writerow({'name': par, 'LaTeX': p[par]})

                # write out parameter definitions for derived model parameters (eg. c0 or sigmav)
                for par in model_params:
                    writer.writerow({'name': par + '*', 'LaTeX': model_params[par]})

        else:

            print ':: GetDist .paramnames file "{f}" already exists: leaving intact'.format(f=getdist_param_file)


        # generate GetDist-compatible chain file from CosmoSIS output, if it does not already exist

        if not os.path.exists(getdist_chain_file):

            print ':: converting CosmoSIS-format chain file "{s}" to GetDist-format chain file "{o}"'.format(s=cosmosis_file, o=getdist_chain_file)

            # note p.keys() must return a list that is ordered in the correct way
            input_columns = list(p.keys()) + list(model_params.keys()) + ['like']
            output_columns = ['weight', 'like'] + list(p.keys()) + list(model_params.keys())

            table = ascii.read(cosmosis_file, Reader=ascii.NoHeader, names=input_columns)

            with open(getdist_chain_file, 'w') as g:

                writer = csv.DictWriter(g, output_columns, delimiter='\t', restval='MISSING')

                for row in table:

                    row_dict = {}

                    # populate with entries for bias model parameters
                    for par in p:
                        row_dict.update({par: row[par]})

                    # populate with entries for derived model parameters (eg. c0 or sigmav)
                    for par in model_params:
                        row_dict.update({par: row[par]})

                    minus_log_lik = -row['like']
                    row_dict.update({'weight': 1, 'like': minus_log_lik})

                    writer.writerow(row_dict)

        else:

            print ':: GetDist-format chain file "{o}" already exists: leaving intact; no conversion of "{s}"'.format(s=cosmosis_file, o=getdist_chain_file)


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

        r = {p: x.parWithName(p).bestfit_sample for p in list(p.keys()) + list(model_params.keys())}

        self.bestfit = r
        self.bestfit_chisquare = self.compute_chisquare(r, mp, glb)


    def get_fit_point(self):

        return self.bestfit


class analyse_maxlike(analyse_core):

    def __init__(self, t, p, root_path, maxlike_file, mp, glb):

        super(analyse_maxlike, self).__init__(t, p)

        # construct hierarchy of plot folders
        plot_folder = os.path.join(root_path, 'plots')

        if not os.path.exists(plot_folder):
            try:
                os.makedirs(plot_folder)
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise


        # READ OUTPUT FILE GENERATED BY MAXLIKE

        # this will consist of a single line giving the best-fit values of all parameters, including derived ones

        # cache OrderedDict of model parameters
        model_params = self.get_model_params()

        # note p.keys() must return a list that is ordered in the correct way
        input_columns = list(p.keys()) + list(model_params.keys()) + ['like']

        table = ascii.read(maxlike_file, Reader=ascii.NoHeader, names=input_columns)

        row = table[0]

        r = {p: row[p] for p in list(p.keys() + list(model_params.keys()))}

        self.bestfit = r
        self.best_chisquare = self.compute_chisquare(r, mp, glb)


    def get_fit_point(self):

        return self.bestfit



def write_summary(analysis_list, root_path, out_file, mp, glb):

    # filenames for GetDist chain-like output
    getdist_param_file = os.path.join(root_path, 'plots', out_file + '.paramnames')
    getdist_chain_file = os.path.join(root_path, 'plots', out_file + '.txt')

    # parameter lists from all realizations should be the same
    params, model_params = __get_parameter_lists(analysis_list)

    with open(getdist_param_file, 'w') as f:

        print ':: generating GetDist-format chain summary .paramnames "{f}"'.format(f=getdist_param_file)

        writer = csv.DictWriter(f, ['name', 'LaTeX'], delimiter='\t', restval='MISSING')

        for p in params:
            writer.writerow({'name': p, 'LaTeX': params[p]})

        for p in model_params:
            writer.writerow({'name': p + '*', 'LaTeX': model_params[p]})

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

        columns = ['weight', 'like'] + list(params.keys()) + list(model_params.keys()) + \
                  ['0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25', '0.27', '0.29']

        writer = csv.DictWriter(f, columns, delimiter='\t', restval='MISSING')

        for real in analysis_list:

            rlz = analysis_list[real]

            row = rlz.get_fit_point()
            tools = rlz.get_tools()

            param_dict = mp(row)
            coeffs = tools.make_coeff_dict(param_dict)
            P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, row, glb(row))

            deviations = tools.compute_chisq_variation(P0, P2, P4)

            row.update(deviations)
            row.update({'weight': 1, 'like': 1})

            writer.writerow(row)


def __get_parameter_lists(analysis_list):

    params = None
    model_params = None

    for r in analysis_list:

        p = analysis_list[r].get_bias_params()
        p_model = analysis_list[r].get_model_params()

        if params is not None:
            if params != p:
                raise RuntimeError
        else:
            params = p

        if model_params is not None:
            if model_params != p_model:
                raise RuntimeError
        else:
            model_params = p_model

    return params, model_params


def write_Pell(list, root_path, out_file, mp, glb):

    for real in list:

        rlz = list[real]

        p = os.path.join(root_path, 'plots', out_file + '_' + real + '_Pell.csv')

        bestfit = rlz.get_fit_point()
        tools = rlz.get_tools()

        params = mp(bestfit)
        coeffs = tools.make_coeff_dict(params)
        P0, P2, P4 = tools.theory.build_theory_P_ell(coeffs, bestfit, glb(bestfit))

        tools.make_summary_plot(p, P0, P2, P4)
