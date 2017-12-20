import numpy as np
from scipy.special import erf

import copy
from collections import OrderedDict

from scipy import optimize

from products import database
from config import settings


# services class for generic exponential (Gaussian) FoG model
class theory(object):

    def __init__(self, my_config, k_sample, col_name):

        # import theory data products in self.data
        self.data = database(my_config, k_sample, col_name)

        self.parameters = OrderedDict([('sigmav', '\sigma_v')])


    def __ell_0(self, k, f, sigma):

        if abs(sigma) > settings.LowSigmaSeries:

            mu0 = (np.sqrt(np.pi)*erf(f*k*sigma))/(2.*f*k*sigma)

            mu2 = ((-2*f*k*sigma)/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + np.sqrt(np.pi)*erf(f*k*sigma))/(4.*np.power(f,3)*np.power(k,3)*np.power(sigma,3))

            mu4 = ((-2*f*k*sigma*(3 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 3*np.sqrt(np.pi)*erf(f*k*sigma))/(8.*np.power(f,5)*np.power(k,5)*np.power(sigma,5))

            mu6 = ((-2*f*k*sigma*(15 + 10*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 15*np.sqrt(np.pi)*erf(f*k*sigma))/(16.*np.power(f,7)*np.power(k,7)*np.power(sigma,7))

            mu8 = ((-2*f*k*sigma*(105 + 70*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 28*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 8*np.power(f,6)*np.power(k,6)*np.power(sigma,6)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 105*np.sqrt(np.pi)*erf(f*k*sigma))/(32.*np.power(f,9)*np.power(k,9)*np.power(sigma,9))

        else:

            mu0 = 1 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/3. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/10.

            mu2 = 0.3333333333333333 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/5. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/14.

            mu4 = 0.2 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/7. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/18.

            mu6 = 0.14285714285714285 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/9. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/22.

            mu8 = 0.1111111111111111 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/11. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/26.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_2(self, k, f, sigma):

        if abs(sigma) > settings.LowSigmaSeries:

            mu0 = (5*((-6*f*k*sigma)/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + np.sqrt(np.pi)*(3 - 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))*erf(f*k*sigma)))/(8.*np.power(f,3)*np.power(k,3)*np.power(sigma,3))

            mu2 = (5*((-2*f*k*sigma*(9 + 4*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + np.sqrt(np.pi)*(9 - 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))*erf(f*k*sigma)))/(16.*np.power(f,5)*np.power(k,5)*np.power(sigma,5))

            mu4 = (5*((-2*f*k*sigma*(45 + 24*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 8*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 3*np.sqrt(np.pi)*(15 - 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))*erf(f*k*sigma)))/(32.*np.power(f,7)*np.power(k,7)*np.power(sigma,7))

            mu6 = (5*((-2*f*k*sigma*(315 + 180*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 64*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 16*np.power(f,6)*np.power(k,6)*np.power(sigma,6)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 15*np.sqrt(np.pi)*(21 - 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))*erf(f*k*sigma)))/(64.*np.power(f,9)*np.power(k,9)*np.power(sigma,9))

            mu8 = (-5*(5670*f*k*sigma + 16*np.power(f,3)*np.power(k,3)*np.power(sigma,3)*(210 + 77*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 20*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 4*np.power(f,6)*np.power(k,6)*np.power(sigma,6)) + 105*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.sqrt(np.pi)*(-27 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))*erf(f*k*sigma)))/(128.*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.power(f,11)*np.power(k,11)*np.power(sigma,11))

        else:

            mu0 = (-2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/3. + (2*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/7.

            mu2 = 0.6666666666666666 - (4*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/7. + (5*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/21.

            mu4 = 0.5714285714285714 - (10*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/21. + (20*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/99.

            mu6 = 0.47619047619047616 - (40*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/99. + (25*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu8 = 0.40404040404040403 - (50*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (2*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/13.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_4(self, k, f, sigma):

        if abs(sigma) > settings.LowSigmaSeries:

            mu0 = (9*((-10*f*k*sigma*(21 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 3*np.sqrt(np.pi)*(35 - 20*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))*erf(f*k*sigma)))/(64.*np.power(f,5)*np.power(k,5)*np.power(sigma,5))

            mu2 = (9*((-2*f*k*sigma*(525 + 170*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 32*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 3*np.sqrt(np.pi)*(175 - 60*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))*erf(f*k*sigma)))/(128.*np.power(f,7)*np.power(k,7)*np.power(sigma,7))

            mu4 = (9*((-2*f*k*sigma*(3675 + 1550*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 416*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 64*np.power(f,6)*np.power(k,6)*np.power(sigma,6)))/np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2)) + 3*np.sqrt(np.pi)*(1225 + 12*np.power(f,2)*np.power(k,2)*np.power(sigma,2)*(-5 + f*k*sigma)*(5 + f*k*sigma))*erf(f*k*sigma)))/(256.*np.power(f,9)*np.power(k,9)*np.power(sigma,9))

            mu6 = (9*(-2*f*k*sigma*(33075 + 15750*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 64*np.power(f,4)*np.power(k,4)*np.power(sigma,4)*(75 + 15*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 2*np.power(f,4)*np.power(k,4)*np.power(sigma,4))) + 45*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.sqrt(np.pi)*(735 - 140*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))*erf(f*k*sigma)))/(512.*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.power(f,11)*np.power(k,11)*np.power(sigma,11))

            mu8 = (9*(-2*f*k*sigma*(363825 + 185850*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 64*np.power(f,4)*np.power(k,4)*np.power(sigma,4)*(945 + 210*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 34*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 4*np.power(f,6)*np.power(k,6)*np.power(sigma,6))) + 315*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.sqrt(np.pi)*(1155 - 180*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))*erf(f*k*sigma)))/(1024.*np.exp(np.power(f,2)*np.power(k,2)*np.power(sigma,2))*np.power(f,13)*np.power(k,13)*np.power(sigma,13))

        else:

            mu0 = (4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/35.

            mu2 = (-8*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/35. + (12*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/77.

            mu4 = 0.22857142857142856 - (24*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/77. + (24*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu6 = 0.3116883116883117 - (48*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (24*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu8 = 0.3116883116883117 - (48*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (24*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __dotp(self, a, b):

        zip = [ a[i]*b[i] for i in xrange(len(a))]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a

        return P


    def build_theory_P_ell(self, coeffs, values, blinear):

        f = self.data.f
        sigmav = values['sigmav']
        
        ks = self.data.k_sample.WiggleZ_conv_ks
        
        ell0_coeff = self.__ell_0(ks, f, sigmav)
        ell2_coeff = self.__ell_2(ks, f, sigmav)
        ell4_coeff = self.__ell_4(ks, f, sigmav)

        zip0 = [self.__dotp(ell0_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]
        zip2 = [self.__dotp(ell2_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]
        zip4 = [self.__dotp(ell4_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]

        P0 = zip0[0].copy()
        for a in zip0[1:]:
            P0 += a  # in-place addition is fastest

        P2 = zip2[0].copy()
        for a in zip2[1:]:
            P2 += a  # in-place addition is fastest

        P4 = zip4[0].copy()
        for a in zip4[1:]:
            P4 += a  # in-place addition is fastest

        return P0, P2, P4


    def compute_model_parameters(self, coeffs, blinear, LikelihoodAgent):
        
        initial_sigmav = np.array([0.0])
        res = optimize.minimize(self.__sigmav_fit, initial_sigmav, method='Powell',
                                args=(coeffs, blinear, LikelihoodAgent),
                                options={'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 50000, 'maxfev': 50000})

        if not res.success:
            raise RuntimeError(res.message)

        return {'sigmav': abs(np.asscalar(res.x))}

    
    def __sigmav_fit(self, x, coeffs, blinear, LikelihoodAgent):

        values = {'sigmav': np.asscalar(x)}

        P0, P2, P4 = self.build_theory_P_ell(coeffs, values, blinear)
        lik = LikelihoodAgent.compute_likelihood(P0, P2, P4, type='ren')

        return -lik

