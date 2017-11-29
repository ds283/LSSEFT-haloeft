import analyse as asy
import coevo_params as coevo

r01 = asy.analyse_emcee(1, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r01.txt', 'coevo_r01',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r02 = asy.analyse_emcee(2, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r02.txt', 'coevo_r02',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r03 = asy.analyse_emcee(3, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r03.txt', 'coevo_r03',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r04 = asy.analyse_emcee(4, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r04.txt', 'coevo_r04',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r05 = asy.analyse_emcee(5, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r05.txt', 'coevo_r05',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r06 = asy.analyse_emcee(6, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r06.txt', 'coevo_r06',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r07 = asy.analyse_emcee(7, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r07.txt', 'coevo_r07',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r08 = asy.analyse_emcee(8, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r08.txt', 'coevo_r08',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r09 = asy.analyse_emcee(9, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r09.txt', 'coevo_r09',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r10 = asy.analyse_emcee(10, {'b1': 'b_1', 'b2': 'b_2'}, 'output_coevo_r10.txt', 'coevo_r10',
                        ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])

list = {'r01': r01, 'r02': r02, 'r03': r03, 'r04': r04, 'r05': r05, 'r06': r06, 'r07': r07, 'r08': r08, 'r09': r09, 'r10': r10}

asy.write_summary(list, coevo.make_params, coevo.get_linear_bias, 'coevo_ensemble')
asy.write_Pell(list, coevo.make_params, coevo.get_linear_bias, 'coevo_ensemble')
