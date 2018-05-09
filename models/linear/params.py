from collections import OrderedDict

param_dict = OrderedDict([('b1', 'b_1')])


def make_params(plist):

    b1 = plist['b1']

    # build coefficient dictionary
    params = {'b1_1': b1, 'b1_2': b1, 'b1_3': b1, 'b2_2': 0.0, 'b2_3': 0.0, 'b3': 0.0,
              'bG2_2': 0.0, 'bG2_3': 0.0, 'bdG2': 0.0, 'bGamma3': 0.0}

    return params


def get_linear_bias(plist):

    return plist['b1']
