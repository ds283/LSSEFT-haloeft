import analyse as asy
import full_params as fp

r01 = asy.analyse_emcee(1, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r01.txt',
                        'plots/full_r01.txt')
r02 = asy.analyse_emcee(2, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r02.txt',
                        'plots/full_r02.txt')
r03 = asy.analyse_emcee(3, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r03.txt',
                        'plots/full_r03.txt')
r04 = asy.analyse_emcee(4, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r04.txt',
                        'plots/full_r04.txt')
r05 = asy.analyse_emcee(5, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r05.txt',
                        'plots/full_r05.txt')
r06 = asy.analyse_emcee(6, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r06.txt',
                        'plots/full_r06.txt')
r07 = asy.analyse_emcee(7, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r07.txt',
                        'plots/full_r07.txt')
r08 = asy.analyse_emcee(8, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r08.txt',
                        'plots/full_r08.txt')
r09 = asy.analyse_emcee(9, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                            'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r09.txt',
                        'plots/full_r09.txt')
r10 = asy.analyse_emcee(10, {'b1_1': 'b_1^{(1)}', 'b1_2': 'b_1^{(2)}', 'b1_3': 'b_1^{(3)}', 'b2_2': 'b_2^{(2)}',
                             'bG2_2': 'b_{G_2}^{(2)}', 'bG2_3': 'b_{G_2}^{(3)}'}, 'output/output_full_r10.txt',
                        'plots/full_r10.txt')

list = {'r01': r01, 'r02': r02, 'r03': r03, 'r04': r04, 'r05': r05, 'r06': r06, 'r07': r07, 'r08': r08, 'r09': r09, 'r10': r10}

asy.write_summary(list, fp.make_params, fp.get_linear_bias, 'full_ensemble')
asy.write_Pell(list, fp.make_params, fp.get_linear_bias, 'full_ensemble')
