import numpy as np

import copy
from collections import OrderedDict

from scipy import optimize

from products import database


# EFT services class
class theory(object):

    def __init__(self, my_config, k_sample):

        # import theory data products into self.data
        self.data = database(my_config, k_sample)

        self.parameters = OrderedDict([('c0', 'c_0'), ('c2', 'c_2'), ('c4', 'c_4'),
                                       ('d1', 'd_1'), ('d2', 'd_2'), ('d3', 'd_3')])


    def build_theory_P_ell(self, coeffs, values, blinear):

        # construct P_ell, including mixing counterterms and stochastic counterterms

        # populate coefficient dictionary with entries mapped from EFT counterterms
        cs = copy.deepcopy(coeffs)
        self.__add_params_to_coeff_dict(cs, values, blinear)

        # this implementation is a bit cryptic but has been tuned for speed -- this is the rate-limiting
        # step in an MCMC analysis
        zip = [cs[key] * data for key, data in self.data.payload.iteritems()]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a  # in-place addition is fastest

        return P[0], P[1], P[2]


    def compute_model_parameters(self, coeffs, blinear, LikelihoodAgent):

        cs = copy.deepcopy(coeffs)

        # optimize fit for counterterms over the renormalization region
        # coeffs is updated with the values of the counterterms
        initial_counterterms = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        res = optimize.minimize(self.__EFT_fit, initial_counterterms, method='Powell',
                                args=(cs, blinear, LikelihoodAgent),
                                options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 50000, 'maxfev': 50000})

        if not res.success:
            raise RuntimeError(res.message)

        return {'c0': res.x[0],
                'c2': res.x[1],
                'c4': res.x[2],
                'd1': res.x[3],
                'd2': res.x[4],
                'd3': res.x[5]}


    def __EFT_fit(self, x, coeffs, blinear, LikelihoodAgent):

        c0, c2, c4, d1, d2, d3 = x

        values = {'c0': c0, 'c2': c2, 'c4': c4, 'd1': d1, 'd2': d2, 'd3': d3}

        P0, P2, P4 = self.build_theory_P_ell(coeffs, values, blinear)
        lik = LikelihoodAgent.compute_likelihood(P0, P2, P4, type='ren')

        return -lik


    def __add_params_to_coeff_dict(self, coeffs, values, blinear):

        fb = self.data.f / blinear

        c0 = values['c0']
        c2 = values['c2']
        c4 = values['c4']
        d1 = values['d1']
        d2 = values['d2']
        d3 = values['d3']

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        # update coefficient dictionary with counterterm values
        coeffs.update({'c0': c0 * blinear,
                       'c2': c2 * blinear,
                       'c4': c4 * blinear,
                       'c6': c6_ope * blinear,
                       'd1': d1,
                       'd2': d2,
                       'd3': d3})
