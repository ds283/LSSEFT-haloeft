import numpy as np

from scipy import optimize


class cosmosis_pipeline(object):

    def __init__(self, my_config, my_name, data, theory):

        self.__mod_name = my_name

        self.data = data
        self.theory = theory


    def make_coeff_dict(self, params):

        # extract parameters
        b1_1 = params['b1_1']
        b1_2 = params['b1_2']
        b1_3 = params['b1_3']
        b2_2 = params['b2_2']
        b2_3 = params['b2_3']
        b3 = params['b3']
        bG2_2 = params['bG2_2']
        bG2_3 = params['bG2_3']
        bdG2 = params['bdG2']
        bGamma3 = params['bGamma3']

        # prepare coefficient dictionary
        coeffs = {'nobias': 1.0,
                  'b1_1': b1_1,
                  'b1_2': b1_2,
                  'b1_3': b1_3,
                  'b2_2': b2_2,
                  'b2_3': b2_3,
                  'b3': b3,
                  'bG2_2': bG2_2,
                  'bG2_3': bG2_3,
                  'bdG2': bdG2,
                  'bGamma3': bGamma3,
                  'b1_1_b1_1': b1_1 * b1_1,
                  'b1_2_b1_2': b1_2 * b1_2,
                  'b1_1_b1_2': b1_1 * b1_2,
                  'b1_1_b1_3': b1_1 * b1_3,
                  'b1_1_b2_2': b1_1 * b2_2,
                  'b1_1_b2_3': b1_1 * b2_3,
                  'b1_2_b2_2': b1_2 * b2_2,
                  'b2_2_b2_2': b2_2 * b2_2,
                  'b1_1_b3': b1_1 * b3,
                  'b1_1_bG2_2': b1_1 * bG2_2,
                  'b1_1_bG2_3': b1_1 * bG2_3,
                  'b1_2_bG2_2': b1_2 * bG2_2,
                  'bG2_2_bG2_2': bG2_2 * bG2_2,
                  'b2_2_bG2_2': b2_2 * bG2_2,
                  'b1_1_bdG2': b1_1 * bdG2,
                  'b1_1_bGamma3': b1_1 * bGamma3}

        return coeffs


    def add_counterterms_to_coeff_dict(self, coeffs, values, blinear):

        fb = self.theory.f / blinear

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


    def compute(self, block, params, blinear, likes):

        coeffs = self.make_coeff_dict(params)

        # compute EFT counterterms
        # note this will update coeffs with the values of counterterms during the optimization process,
        # but these will be overwritten with the bestfit values immediately below
        values = self.__compute_EFT_counterterms(coeffs, blinear)
        self.add_counterterms_to_coeff_dict(coeffs, values, blinear)

        # build theoretical P_ell with these counterterms
        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        # sum likelihood over all regions and store back into the datablock
        lik = sum([self.compute_likelihood(region, P0, P2, P4, 'fit') for region in self.data.regions])
        block[likes, 'HALOEFT_LIKE'] = lik


        # store derived EFT parameters for output with the rest of the chain
        block['counterterms', 'c0'] = values['c0']
        block['counterterms', 'c2'] = values['c2']
        block['counterterms', 'c4'] = values['c4']

        block['counterterms', 'd1'] = values['d1']
        block['counterterms', 'd2'] = values['d2']
        block['counterterms', 'd3'] = values['d3']

        return 0


    def build_theory_P_ell(self, coeffs):

        # construct P_ell, including mixing counterterms and stochastic counterterms

        # this implementation is a bit cryptic but has been tuned for speed -- this is the rate-limiting
        # step in an MCMC analysis
        zip = [ coeffs[key] * data for key, data in self.theory.payload.iteritems() ]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a              # in-place addition is fastest

        return P[0], P[1], P[2]


    def compute_likelihood(self, region, P0, P2, P4, type='fit'):

        if type is 'fit':
            mask = self.data.k_sample.conv_fit_mask
            means = self.data.fit_means[region]
            inv_cov = self.data.fit_inv_covs[region]
        else:
            mask = self.data.k_sample.conv_ren_mask
            means = self.data.ren_means[region]
            inv_cov = self.data.ren_inv_covs[region]

        conv = self.data.convs[region]

        # first, convolve theory vector with WizCOLA convolution matrix for this region
        P = np.concatenate( (P0, P2, P4) )
        Pconv = np.dot(conv, P)

        # strip out components that can be compared with measurement
        Ptheory = Pconv[mask]

        # compute chi^2
        Delta = Ptheory - means
        return -np.dot(np.dot(Delta, inv_cov), Delta) / 2.0


    def __compute_EFT_counterterms(self, coeffs, blinear):

        fb = self.theory.f / blinear

        # optimize fit for counterterms over the renormalization region
        # coeffs is updated with the values of the counterterms
        initial_counterterms = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        res = optimize.minimize(self.__compute_renormalization_fit, initial_counterterms, method='Powell',
                                args=(coeffs, blinear, fb), options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 50000, 'maxfev': 50000})

        if not res.success:
            raise RuntimeError(res.message)

        return {'c0': res.x[0],
                'c2': res.x[1],
                'c4': res.x[2],
                'd1': res.x[3],
                'd2': res.x[4],
                'd3': res.x[5]}


    def __compute_renormalization_fit(self, x, coeffs, blinear, fb):

        c0, c2, c4, d1, d2, d3 = x

        # OPE constrains coefficient of mu^6 counterterm
        c6_ope = fb*c4 - fb*fb*c2 + fb*fb*fb*c0

        coeffs.update({'c0': c0 * blinear,
                       'c2': c2 * blinear,
                       'c4': c4 * blinear,
                       'c6': c6_ope * blinear,
                       'd1': d1,
                       'd2': d2,
                       'd3': d3})

        P0, P2, P4 = self.build_theory_P_ell(coeffs)

        lik = sum([self.compute_likelihood(region, P0, P2, P4, 'ren') for region in self.data.regions])

        return -lik


    def cleanup(self):

        return 0
