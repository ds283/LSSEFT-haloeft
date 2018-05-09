from collections import OrderedDict
import multiprocessing as mp
import traceback
import argparse
import os
import imp

import WizCOLA

from RedshiftModels import EFT
from RedshiftModels import Kaiser
from RedshiftModels import OneLoop

import haloeft as heft

from analysis import analyse as asy
from analysis.config import make_config_block


realizations = ['ensemble', 'r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10']

realization_numbers = {'ensemble': 0, 'r01': 1, 'r02': 2, 'r03': 3, 'r04': 4, 'r05': 5,
                       'r06': 6, 'r07': 7, 'r08': 8, 'r09': 9, 'r10': 10}

inputs = {'ensemble': 'ensemble.txt',
          'r01': 'r01.txt', 'r02': 'r02.txt', 'r03': 'r03.txt',
          'r04': 'r04.txt', 'r05': 'r05.txt', 'r06': 'r06.txt',
          'r07': 'r07.txt', 'r08': 'r08.txt', 'r09': 'r09.txt',
          'r10': 'r10.txt'}

# outputs should be GetDist-format roots, which *don't* include the .txt extension
outputs = {'ensemble': 'ensemble',
           'r01': 'r01', 'r02': 'r02', 'r03': 'r03', 'r04': 'r04',
           'r05': 'r05', 'r06': 'r06', 'r07': 'r07', 'r08': 'r08',
           'r09': 'r09', 'r10': 'r10'}

mixing_params = None
stochastic_params = None


def f((realization_tag, file_name, root_path, params_module)):

    try:

        if 'EFT' in root_path:
            rsd_model = 'EFT'
        elif 'KaiserTree' in root_path:
            rsd_model = 'KaiserTree'
        elif 'KaiserHalofit' in root_path:
            rsd_model = 'KaiserHalofit'
        elif 'OneLoop' in root_path:
            rsd_model = 'OneLoop'
        else:
            print "Cannot deduce RSD model used by '{p}'".format(p=root_path)
            raise RuntimeError

        config = make_config_block(realization_numbers[realization_tag], rsd_model)

        ks = WizCOLA.ksamples(config)
        data = WizCOLA.products(config, ks)

        if rsd_model is 'EFT':
            theory = EFT.theory(config, ks)
        elif rsd_model is 'KaiserTree' or rsd_model is 'KaiserHalofit':
            theory = Kaiser.theory(config, ks)
        elif rsd_model is 'OneLoop':
            theory = OneLoop.theory(config, ks)
        else:
            print "Unknown RSD model label '{p}'".format(p=rsd_model)
            raise RuntimeError

        t = heft.tools(root_path, data, theory)

        ptools = imp.load_source("params", params_module)

        obj = asy.analyse_cosmosis(t, ptools.param_dict, root_path, file_name,
                                   outputs[realization_tag], outputs['ensemble'],
                                   ptools.make_params, ptools.get_linear_bias,
                                   mixing_params, stochastic_params)

    except Exception as e:

        print 'Caught exception in worker thread'

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        raise e

    return realization_tag, obj


if __name__ == '__main__':

    # build command-line parser to extract pathnames from the command line
    parser = argparse.ArgumentParser(description='Perform post-processing analysis of halo fits')
    parser.add_argument('dirnames', metavar='path', type=str, nargs='*',
                        help='path name containing CosmoSIS output products')

    args = parser.parse_args()

    # loop through each pathname provided
    for path in args.dirnames:

        # check whether path exists
        if not os.path.exists(path):

            print "ignoring non-existent path '{p}'".format(p=path)
            continue

        # check with params module is present
        params_module = os.path.join(path, '..', 'params.py')
        if not os.path.exists(params_module):

            print "ignoring path '{p}'; cannot find parameters module params.py".format(p=path)
            continue

        print "-- analysing path '{p}'".format(p=path)

        # check which region files exist
        output_path = os.path.join(path, 'output')

        realization_files = []

        for real in realizations:

            file = os.path.join(output_path, inputs[real])

            if os.path.exists(file):

                realization_files.append((real, file, path, params_module))
                print "** added realization output '{p}' to analysis list".format(p=file)

        if len(realization_files) > 0:

            # list will hold the analysis objects for each file
            obj_list = OrderedDict()

            # set up a multiprocessing pool to parallelize the analysis step
            p = mp.Pool()

            # build analysis objects for each file we detected
            for real, obj in p.map(f, realization_files):

                obj_list[real] = obj

            p.close()

            ptools = imp.load_source("params", params_module)

            # write summary CSV file of best-fit parameter values and chi-squares (as a function of k)
            asy.write_summary(obj_list, path, 'bestfit_chisq', ptools.make_params, ptools.get_linear_bias)

            # write more detailed CSV filea containing best-fit power spectra and WizCOLA data points,
            # both for best-fit point and ensemble-average point (if exists)
            asy.write_bestfit_plots(obj_list, path, 'bestfit', 'ensemblefit', ptools.make_params, ptools.get_linear_bias)

        else:

            print '** did not find any realization outputs to analyse'
