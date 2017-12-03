import analyse as asy
from collections import OrderedDict
import linear_params as lin

params = OrderedDict([('b1', 'b_1')])

r01 = asy.analyse_emcee(1, params, 'output_linear_r01.txt', 'linear_r01', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r02 = asy.analyse_emcee(2, params, 'output_linear_r02.txt', 'linear_r02', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r03 = asy.analyse_emcee(3, params, 'output_linear_r03.txt', 'linear_r03', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r04 = asy.analyse_emcee(4, params, 'output_linear_r04.txt', 'linear_r04', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r05 = asy.analyse_emcee(5, params, 'output_linear_r05.txt', 'linear_r05', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r06 = asy.analyse_emcee(6, params, 'output_linear_r06.txt', 'linear_r06', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r07 = asy.analyse_emcee(7, params, 'output_linear_r07.txt', 'linear_r07', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r08 = asy.analyse_emcee(8, params, 'output_linear_r08.txt', 'linear_r08', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r09 = asy.analyse_emcee(9, params, 'output_linear_r09.txt', 'linear_r09', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])
r10 = asy.analyse_emcee(10, params, 'output_linear_r10.txt', 'linear_r10', lin.make_params, lin.get_linear_bias,
                        ['b1', 'c0', 'c2', 'c4'],
                        ['b1', 'd1', 'd2', 'd3'])

list = OrderedDict(
    [('r01', r01), ('r02', r02), ('r03', r03), ('r04', r04), ('r05', r05), ('r06', r06), ('r07', r07), ('r08', r08),
     ('r09', r09), ('r10', r10)])

asy.write_summary(list, 'linear_ensemble')
asy.write_Pell(list, 'linear_ensemble')
