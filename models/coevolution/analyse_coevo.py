from analysis import analyse as asy
from collections import OrderedDict
import multiprocessing as mp
import traceback
import coevo_params as param_tools

tag = 'coevo'
model_name = 'Coevolution'

params = OrderedDict([('b1', 'b_1'), ('b2', 'b_2')])

regions = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10']

numbers = {'r01': 1, 'r02': 2, 'r03': 3, 'r04': 4, 'r05': 5, 'r06': 6, 'r07': 7, 'r08': 8, 'r09': 9, 'r10': 10}

inputs = {'r01': 'output_'+tag+'_r01.txt', 'r02': 'output_'+tag+'_r02.txt', 'r03': 'output_'+tag+'_r03.txt',
          'r04': 'output_'+tag+'_r04.txt', 'r05': 'output_'+tag+'_r05.txt', 'r06': 'output_'+tag+'_r06.txt',
          'r07': 'output_'+tag+'_r07.txt', 'r08': 'output_'+tag+'_r08.txt', 'r09': 'output_'+tag+'_r09.txt',
          'r10': 'output_'+tag+'_r10.txt'}

outputs = {'r01': tag+'_r01', 'r02': tag+'_r02', 'r03': tag+'_r03', 'r04': tag+'_r04',
           'r05': tag+'_r05', 'r06': tag+'_r06', 'r07': tag+'_r07', 'r08': tag+'_r08',
           'r09': tag+'_r09', 'r10': tag+'_r10'}

mixing_params = ['b1', 'b2', 'c0', 'c2', 'c4']
stochastic_params = ['b1', 'b2', 'd1', 'd2', 'd3']


def f(tag):

    try:

        obj = asy.analyse_CosmoSIS(numbers[tag], model_name, params, inputs[tag], outputs[tag],
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

    asy.write_summary(list, tag+'_ensemble')
    asy.write_Pell(list, tag+'_ensemble')
