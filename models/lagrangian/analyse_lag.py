import analyse as asy
from collections import OrderedDict
import lag_params as lag

params = OrderedDict([('b1', 'b_1'), ('b2', 'b_2'), ('bs2', 'b_{s^2}'), ('b3nl', r'b_{3\mathrm{nl}}')])

r01 = asy.analyse_emcee(1, params, 'output_lag_r01.txt', 'lag_r01', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r02 = asy.analyse_emcee(2, params, 'output_lag_r02.txt', 'lag_r02', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r03 = asy.analyse_emcee(3, params, 'output_lag_r03.txt', 'lag_r03', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r04 = asy.analyse_emcee(4, params, 'output_lag_r04.txt', 'lag_r04', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r05 = asy.analyse_emcee(5, params, 'output_lag_r05.txt', 'lag_r05', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r06 = asy.analyse_emcee(6, params, 'output_lag_r06.txt', 'lag_r06', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r07 = asy.analyse_emcee(7, params, 'output_lag_r07.txt', 'lag_r07', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r08 = asy.analyse_emcee(8, params, 'output_lag_r08.txt', 'lag_r08', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r09 = asy.analyse_emcee(9, params, 'output_lag_r09.txt', 'lag_r09', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r10 = asy.analyse_emcee(10, params, 'output_lag_r10.txt', 'lag_r10', lag.make_params, lag.get_linear_bias,
                        ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])

list = OrderedDict(
    [('r01', r01), ('r02', r02), ('r03', r03), ('r04', r04), ('r05', r05), ('r06', r06), ('r07', r07), ('r08', r08),
     ('r09', r09), ('r10', r10)])

asy.write_summary(list, 'lag_ensemble')
asy.write_Pell(list, 'lag_ensemble')
