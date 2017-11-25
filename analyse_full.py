import analyse as asy
from models.full import full_params as fp

r01 = asy.analyse(1, 'output/output_full_r01.txt', 'plots/full_r01.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r02 = asy.analyse(2, 'output/output_full_r02.txt', 'plots/full_r02.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r03 = asy.analyse(3, 'output/output_full_r03.txt', 'plots/full_r03.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r04 = asy.analyse(4, 'output/output_full_r04.txt', 'plots/full_r04.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r05 = asy.analyse(5, 'output/output_full_r05.txt', 'plots/full_r05.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r06 = asy.analyse(6, 'output/output_full_r06.txt', 'plots/full_r06.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r07 = asy.analyse(7, 'output/output_full_r07.txt', 'plots/full_r07.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r08 = asy.analyse(8, 'output/output_full_r08.txt', 'plots/full_r08.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r09 = asy.analyse(9, 'output/output_full_r09.txt', 'plots/full_r09.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})
r10 = asy.analyse(10, 'output/output_full_r10.txt', 'plots/full_r10.txt', {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}', 'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'})

list = {'r01': r01, 'r02': r02, 'r03': r03, 'r04': r04, 'r05': r05, 'r06': r06, 'r07': r07, 'r08': r08, 'r09': r09, 'r10': r10}

asy.write_summary(list, fp.make_params, fp.get_linear_bias, 'full_ensemble')
asy.write_Pell(list, fp.make_params, fp.get_linear_bias, 'full_ensemble')
