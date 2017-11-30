import analyse as asy
from collections import OrderedDict
import lag_params as lag

r01 = asy.analyse_emcee(1, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r01.txt', 'lag_r01', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r02 = asy.analyse_emcee(2, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r02.txt', 'lag_r02', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r03 = asy.analyse_emcee(3, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r03.txt', 'lag_r03', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r04 = asy.analyse_emcee(4, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r04.txt', 'lag_r04', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r05 = asy.analyse_emcee(5, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r05.txt', 'lag_r05', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r06 = asy.analyse_emcee(6, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r06.txt', 'lag_r06', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r07 = asy.analyse_emcee(7, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r07.txt', 'lag_r07', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r08 = asy.analyse_emcee(8, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r08.txt', 'lag_r08', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r09 = asy.analyse_emcee(9, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r09.txt', 'lag_r09', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])
r10 = asy.analyse_emcee(10, {'b1': 'b_1', 'b2': 'b_2', 'bs2': 'b_{s^2}', 'b3nl': r'b_{3\text{nl}}'},
                        'output_lag_r10.txt', 'lag_r10', ['b1', 'b2', 'bs2', 'b3nl', 'c0', 'c2', 'c4'],
                        ['b1', 'b2', 'bs2', 'b3nl', 'd1', 'd2', 'd3'])

list = OrderedDict(
    [('r01', r01), ('r02', r02), ('r03', r03), ('r04', r04), ('r05', r05), ('r06', r06), ('r07', r07), ('r08', r08),
     ('r09', r09), ('r10', r10)])

asy.write_summary(list, lag.make_params, lag.get_linear_bias, 'lag_ensemble')
asy.write_Pell(list, lag.make_params, lag.get_linear_bias, 'lag_ensemble')
