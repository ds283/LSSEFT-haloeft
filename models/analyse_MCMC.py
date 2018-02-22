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


regions = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10']

numbers = {'r01': 1, 'r02': 2, 'r03': 3, 'r04': 4, 'r05': 5, 'r06': 6, 'r07': 7, 'r08': 8, 'r09': 9, 'r10': 10}

inputs = {'r01': 'r01.txt', 'r02': 'r02.txt', 'r03': 'r03.txt',
          'r04': 'r04.txt', 'r05': 'r05.txt', 'r06': 'r06.txt',
          'r07': 'r07.txt', 'r08': 'r08.txt', 'r09': 'r09.txt',
          'r10': 'r10.txt'}

outputs = {'r01': 'r01', 'r02': 'r02', 'r03': 'r03', 'r04': 'r04',
           'r05': 'r05', 'r06': 'r06', 'r07': 'r07', 'r08': 'r08',
           'r09': 'r09', 'r10': 'r10'}

mixing_params = None
stochastic_params = None


def f((region_tag, file_name, root_path, params_module)):

    try:

        if 'Halofit' in root_path:
            config = make_config_block(numbers[region_tag], True)
        else:
            config = make_config_block(numbers[region_tag], False)

        ks = WizCOLA.ksamples(config)
        data = WizCOLA.products(config, ks)

        if 'EFT' in root_path:
            theory = EFT.theory(config, ks)
        elif 'Kaiser' in root_path:
            theory = Kaiser.theory(config, ks)
        elif 'OneLoop' in root_path or 'OneLoop' in root_path:
            theory = OneLoop.theory(config, ks)
        else:
            print "Cannot deduce RSD model used by '{p}'".format(p=root_path)
            raise RuntimeError

        t = heft.tools(root_path, data, theory)

        ptools = imp.load_source("params", params_module)

        obj = asy.analyse_cosmosis(t, ptools.param_dict, root_path, file_name, outputs[region_tag],
                                   ptools.make_params, ptools.get_linear_bias,
                                   mixing_params, stochastic_params)

    except Exception as e:

        print 'Caught exception in worker thread'

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        raise e

    return region_tag, obj


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
        region_files = []

        for r in regions:

            file = os.path.join(output_path, inputs[r])

            if os.path.exists(file):

                region_files.append((r, file, path, params_module))
                print "** added region output '{p}' to analysis list".format(p=file)

        if len(region_files) > 0:

            # list will hold the analysis objects for each file
            obj_list = OrderedDict()

            # set up a multiprocessing pool to parallelize the analysis step
            p = mp.Pool()

            # build analysis objects for each file we detected
            for r, obj in p.map(f, region_files):

                obj_list[r] = obj

            p.close()

            ptools = imp.load_source("params", params_module)

            asy.write_summary(obj_list, path, 'ensemble', ptools.make_params, ptools.get_linear_bias)
            asy.write_Pell(obj_list, path, 'ensemble', ptools.make_params, ptools.get_linear_bias)

        else:

            print '** did not find any region outputs to analyse'
