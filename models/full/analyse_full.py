import analyse as asy
from collections import OrderedDict
import full_params as fp

params = OrderedDict([('b1_1', 'b_1^{(1)}'), ('b1_2', 'b_1^{(2)}'), ('b1_3', 'b_1^{(3)}'),
                      ('b2_2', 'b_2^{(2)}'), ('bG2_2', 'b_{G_2}^{(2)}'), ('bG2_3', 'b_{G_2}^{(3)}')])

r01 = asy.analyse_emcee(1, params, 'output/output_full_r01.txt', 'plots/full_r01.txt', fp.make_params, fp.get_linear_bias)
r02 = asy.analyse_emcee(2, params, 'output/output_full_r02.txt', 'plots/full_r02.txt', fp.make_params, fp.get_linear_bias)
r03 = asy.analyse_emcee(3, params, 'output/output_full_r03.txt', 'plots/full_r03.txt', fp.make_params, fp.get_linear_bias)
r04 = asy.analyse_emcee(4, params, 'output/output_full_r04.txt', 'plots/full_r04.txt', fp.make_params, fp.get_linear_bias)
r05 = asy.analyse_emcee(5, params, 'output/output_full_r05.txt', 'plots/full_r05.txt', fp.make_params, fp.get_linear_bias)
r06 = asy.analyse_emcee(6, params, 'output/output_full_r06.txt', 'plots/full_r06.txt', fp.make_params, fp.get_linear_bias)
r07 = asy.analyse_emcee(7, params, 'output/output_full_r07.txt', 'plots/full_r07.txt', fp.make_params, fp.get_linear_bias)
r08 = asy.analyse_emcee(8, params, 'output/output_full_r08.txt', 'plots/full_r08.txt', fp.make_params, fp.get_linear_bias)
r09 = asy.analyse_emcee(9, params, 'output/output_full_r09.txt', 'plots/full_r09.txt', fp.make_params, fp.get_linear_bias)
r10 = asy.analyse_emcee(10, params, 'output/output_full_r10.txt', 'plots/full_r10.txt', fp.make_params, fp.get_linear_bias)

list = OrderedDict(
    [('r01', r01), ('r02', r02), ('r03', r03), ('r04', r04), ('r05', r05), ('r06', r06), ('r07', r07), ('r08', r08),
     ('r09', r09), ('r10', r10)])

asy.write_summary(list, 'full_ensemble')
asy.write_Pell(list, 'full_ensemble')
