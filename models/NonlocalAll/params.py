from collections import OrderedDict

param_dict = OrderedDict([('b1_1', 'b_1^{(1)}'),
                          ('b1_2', 'b_1^{(2)}'),
                          ('b1_3', 'b_1^{(3)}'),
                          ('b2_2', 'b_2^{(2)}'),
                          ('b2_3', 'b_2^{(3)}'),
                          ('b3', 'b_3'),
                          ('bG2_2', 'b_{G_2}^{(2)}'),
                          ('bG2_3', 'b_{G_2}^{(3)}'),
                          ('bdG2', 'b_{\\delta G_2}'),
                          ('bGamma3', 'b_{\\Gamma_3}')])


def make_params(plist):

    return plist


def get_linear_bias(plist):

    return plist['b1_1']