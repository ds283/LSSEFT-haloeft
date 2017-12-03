import analyse as asy
from collections import OrderedDict
import coevo_params as coevo

params = OrderedDict([('b1', 'b_1'), ('b2', 'b_2')])

r01 = asy.analyse_emcee(1, params, 'output_coevo_r01.txt', 'coevo_r01', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r02 = asy.analyse_emcee(2, params, 'output_coevo_r02.txt', 'coevo_r02', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r03 = asy.analyse_emcee(3, params, 'output_coevo_r03.txt', 'coevo_r03', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r04 = asy.analyse_emcee(4, params, 'output_coevo_r04.txt', 'coevo_r04', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r05 = asy.analyse_emcee(5, params, 'output_coevo_r05.txt', 'coevo_r05', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r06 = asy.analyse_emcee(6, params, 'output_coevo_r06.txt', 'coevo_r06', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r07 = asy.analyse_emcee(7, params, 'output_coevo_r07.txt', 'coevo_r07', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r08 = asy.analyse_emcee(8, params, 'output_coevo_r08.txt', 'coevo_r08', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r09 = asy.analyse_emcee(9, params, 'output_coevo_r09.txt', 'coevo_r09', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r10 = asy.analyse_emcee(10, params, 'output_coevo_r10.txt', 'coevo_r10', coevo.make_params, coevo.get_linear_bias,
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])

list = OrderedDict(
    [('r01', r01), ('r02', r02), ('r03', r03), ('r04', r04), ('r05', r05), ('r06', r06), ('r07', r07), ('r08', r08),
     ('r09', r09), ('r10', r10)])

asy.write_summary(list, 'coevo_ensemble')
asy.write_Pell(list, 'coevo_ensemble')
