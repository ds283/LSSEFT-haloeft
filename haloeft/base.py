import numpy as np
import cpuinfo


class base(object):

    def __init__(self, data, theory):

        self.data = data
        self.theory = theory

        cpu_data = cpuinfo.get_cpu_info()

        print 'LSSEFT-haloeft running on {brand} at {hz_real} (advertised {hz_adv})'.format(
            brand=cpu_data['brand'], hz_real=cpu_data['hz_actual'], hz_adv=cpu_data['hz_advertised'])


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


    def compute_likelihood(self, P0, P2, P4, type='fit'):

        # sum likelihood over all realizations
        lik = sum([self.__compute_region_likelihood(region, P0, P2, P4, type) for region in self.data.regions])

        return lik


    def __compute_region_likelihood(self, region, P0, P2, P4, type='fit'):

        # set up an appropriate set of masks and data items

        if type is 'fit':

            mask = self.data.k_sample.conv_fit_mask
            means = self.data.fit_means[region]
            inv_cov = self.data.fit_inv_covs[region]

        elif type is 'ren':

            mask = self.data.k_sample.conv_ren_mask
            means = self.data.ren_means[region]
            inv_cov = self.data.ren_inv_covs[region]

        else:

            print 'unknown likelihood type'
            raise RuntimeError

        conv = self.data.convs[region]

        # first, convolve theory vector with WizCOLA convolution matrix for this region
        P = np.concatenate( (P0, P2, P4) )
        Pconv = np.dot(conv, P)

        # strip out components that can be compared with measurement
        Ptheory = Pconv[mask]

        # compute chi^2
        Delta = Ptheory - means
        return -np.dot(np.dot(Delta, inv_cov), Delta) / 2.0
