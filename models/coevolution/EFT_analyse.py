from collections import OrderedDict
import multiprocessing as mp
import traceback

import params as param_tools

import WizCOLA
import EFT
import haloeft as heft

from analysis import analyse as asy
from analysis.config import make_config_block


model_name = 'Coevolution'

params = OrderedDict([('b1', 'b_1'), ('b2', 'b_2')])

regions = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10']

numbers = {'r01': 1, 'r02': 2, 'r03': 3, 'r04': 4, 'r05': 5, 'r06': 6, 'r07': 7, 'r08': 8, 'r09': 9, 'r10': 10}

inputs = {'r01': 'EFT_r01.txt', 'r02': 'EFT_r02.txt', 'r03': 'EFT_r03.txt',
          'r04': 'EFT_r04.txt', 'r05': 'EFT_r05.txt', 'r06': 'EFT_r06.txt',
          'r07': 'EFT_r07.txt', 'r08': 'EFT_r08.txt', 'r09': 'EFT_r09.txt',
          'r10': 'EFT_r10.txt'}

outputs = {'r01': 'EFT_r01', 'r02': 'EFT_r02', 'r03': 'EFT_r03', 'r04': 'EFT_r04',
           'r05': 'EFT_r05', 'r06': 'EFT_r06', 'r07': 'EFT_r07', 'r08': 'EFT_r08',
           'r09': 'EFT_r09', 'r10': 'EFT_r10'}

mixing_params = ['b1', 'b2', 'c0', 'c2', 'c4']
stochastic_params = ['b1', 'b2', 'd1', 'd2', 'd3']


def f(tag):

    try:

        config = make_config_block(numbers[tag])
        ks = WizCOLA.ksamples(config)
        data = WizCOLA.products(config, ks)
        theory = EFT.theory(config, ks)

        t = heft.tools(model_name, data, theory)

        obj = asy.analyse_cosmosis(t, params, inputs[tag], outputs[tag],
                                   param_tools.make_params, param_tools.get_linear_bias,
                                   mixing_params, stochastic_params)

    except Exception as e:

        print 'Caught exception in worker thread'

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        raise e

    return obj


if __name__ == '__main__':

    list = OrderedDict()

    p = mp.Pool()

    for n, r in enumerate(p.map(f, regions)):

        label = regions[n]
        list[label] = r

    p.close()

    asy.write_summary(list, 'EFT_ensemble')
    asy.write_Pell(list, 'EFT_ensemble')
