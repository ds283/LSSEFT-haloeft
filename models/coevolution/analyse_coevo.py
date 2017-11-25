import analyse as asy
from models.coevolution import coevo_params as coevo

r01 = asy.analyse(1, 'output_LL_r01.txt', 'LL_r01', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r02 = asy.analyse(2, 'output_LL_r02.txt', 'LL_r02', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r03 = asy.analyse(3, 'output_LL_r03.txt', 'LL_r03', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r04 = asy.analyse(4, 'output_LL_r04.txt', 'LL_r04', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r05 = asy.analyse(5, 'output_LL_r05.txt', 'LL_r05', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r06 = asy.analyse(6, 'output_LL_r06.txt', 'LL_r06', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r07 = asy.analyse(7, 'output_LL_r07.txt', 'LL_r07', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r08 = asy.analyse(8, 'output_LL_r08.txt', 'LL_r08', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r09 = asy.analyse(9, 'output_LL_r09.txt', 'LL_r09', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])
r10 = asy.analyse(10, 'output_LL_r10.txt', 'LL_r10', {'b1': 'b_1', 'b2': 'b_2'}, ['b1', 'b2', 'c0', 'c2', 'c4'], ['b1', 'b2', 'd1', 'd2', 'd3'])

list = {'r01': r01, 'r02': r02, 'r03': r03, 'r04': r04, 'r05': r05, 'r06': r06, 'r07': r07, 'r08': r08, 'r09': r09, 'r10': r10}

asy.write_summary(list, coevo.make_params, coevo.get_linear_bias, 'LL_ensemble')
asy.write_Pell(list, coevo.make_params, coevo.get_linear_bias, 'LL_ensemble')
