from collections import OrderedDict

param_dict = OrderedDict([('b1_1', 'b_1^{(1)}'), ('b1_2', 'b_1^{(2)}'), ('b1_3', 'b_1^{(3)}'),
                          ('b2_2', 'b_2^{(2)}'), ('bG2_2', 'b_{G_2}^{(2)}'), ('bG2_3', 'b_{G_2}^{(3)}')])


def make_params(plist):

    b1_1 = plist['b1_1']
    b1_2 = plist['b1_2']
    b1_3 = plist['b1_3']
    b2_2 = plist['b2_2']
    bG2_2 = plist['bG2_2']
    bG2_3 = plist['bG2_3']

    params = {'b1_1': b1_1, 'b1_2': b1_2, 'b1_3': b1_3, 'b2_2': b2_2, 'b2_3': 0.0, 'b3': 0.0,
              'bG2_2': bG2_2, 'bG2_3': bG2_3, 'bdG2': 0.0, 'bGamma3': 0.0}

    return params


def get_linear_bias(plist):

    return plist['b1_1']
