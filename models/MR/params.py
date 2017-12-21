from collections import OrderedDict

param_dict = OrderedDict([('b1', 'b_1'), ('b2', 'b_2'), ('bs2', 'b_{s^2}'), ('b3nl', r'b_{3\mathrm{nl}}')])


def make_params(plist):

    b1 = plist['b1']
    b2 = plist['b2']
    bs2 = plist['bs2']
    b3nl = plist['b3nl']

    # the local model involves the combination bs2*s2/2 + 105*b3nl*(st + 8*delta*delta*delta/189)/16

    # now, s2 is defined as G2 + 2*delta*delta/3, so if bG2 is defined without a factor of 1/2 and
    # b2 is defined with a factor of 1/2
    bG2 = bs2 / 2.0
    b2 += 2.0 * bs2 / 3.0

    # and st is defined as (the third order part of) Gamma3/2 + delta*delta/3 - theta*theta/3
    # using (3.6) of Assassi et al. arXiv 1402.5916, in the Einstein-de Sitter approximation
    # delta(2) - theta(2) = -(2/7) G2(2)
    # therefore delta*delta/3 - theta*delta/3 = (2/3)[ delta(1)delta(2) - theta(1)theta(2) ]
    # = - (4/21) delta(1) G2(2)
    # so st = Gamma3/2 - (4/21) delta G2 in the EdS approximation (we don't strictly use this,
    # but we suppose this relation is still a reasonable approximation)
    # [note, delta = delta(1) + delta(2) + ... in Assassi et al.; see their Eq. (2.15)]
    bGamma3 = 105.0 * b3nl / 32.0

    # net factor of bdG2 = (105/16)(-4/21) = -5/4
    bdG2 = -5.0 * b3nl / 4.0

    # b3 is defined with a factor of 1/6, so its net factor is 105 * (8/189) * (6/16) = 5/3
    b3 = 5.0 * b3nl / 3.0

    # build coefficient dictionary
    params = {'b1_1': b1, 'b1_2': b1, 'b1_3': b1, 'b2_2': b2, 'b2_3': b2, 'b3': b3,
              'bG2_2': bG2, 'bG2_3': bG2, 'bdG2': bdG2, 'bGamma3': bGamma3}

    return params


def get_linear_bias(plist):

    return plist['b1']
