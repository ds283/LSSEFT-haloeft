import numpy as np


class cosmosis_pipeline(object):

    def __init__(self, my_name, data, theory):

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


    def compute(self, block, params, blinear, likes):

        coeffs = self.make_coeff_dict(params)

        # Most models will depend on extra parameters beyond the bias model, such as EFT counterterms
        # or the velocity dispersion in a Gaussian FoG mode, which we may want to optimize.
        # The theory bundle provides a compute_model_parameters() method for this purpose
        values = self.theory.compute_model_parameters(coeffs, blinear, self)

        # build theoretical P_ell with these bias coefficients and model values
        P0, P2, P4 = self.theory.build_theory_P_ell(coeffs, values, blinear)

        # sum likelihood over all regions and store back into the datablock
        block[likes, 'HALOEFT_LIKE'] = self.compute_likelihood(P0, P2, P4, type='fit')

        # store derived model parameters for output with the rest of the chain
        for p in values:
            block['counterterms', p] = values[p]

        return 0


    def compute_likelihood(self, P0, P2, P4, type='fit'):

        lik = sum([self.__compute_region_likelihood(region, P0, P2, P4, type) for region in self.data.regions])

        return lik


    def __compute_region_likelihood(self, region, P0, P2, P4, type='fit'):

        # set up an appropriate
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


    def cleanup(self):

        return 0
