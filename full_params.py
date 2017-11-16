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


def get_blinear(plist):
    return plist['b1_1']
