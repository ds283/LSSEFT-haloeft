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

            expf = np.exp(-np.power(f,2)*np.power(k,2)*np.power(sigma,2))
            erff = erf(f*k*sigma)

            with np.errstate(over='raise'):

                try:

                    mu0 = expf*(0) + erff*(np.sqrt(np.pi)/(2.*f*k*sigma))

                    mu2 = expf*(-1/(2.*np.power(f,2)*np.power(k,2)*np.power(sigma,2))) + erff*(np.sqrt(np.pi)/(4.*np.power(f,3)*np.power(k,3)*np.power(sigma,3)))

                    mu4 = expf*(-(3 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/(4.*np.power(f,4)*np.power(k,4)*np.power(sigma,4))) + erff*((3*np.sqrt(np.pi))/(8.*np.power(f,5)*np.power(k,5)*np.power(sigma,5)))

                    mu6 = expf*(-(15 + 10*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/(8.*np.power(f,6)*np.power(k,6)*np.power(sigma,6))) + erff*((15*np.sqrt(np.pi))/(16.*np.power(f,7)*np.power(k,7)*np.power(sigma,7)))

                    mu8 = expf*(-(105 + 70*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 28*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 8*np.power(f,6)*np.power(k,6)*np.power(sigma,6))/(16.*np.power(f,8)*np.power(k,8)*np.power(sigma,8))) + erff*((105*np.sqrt(np.pi))/(32.*np.power(f,9)*np.power(k,9)*np.power(sigma,9)))

                except FloatingPointError:

                    print 'k = {k}, f = {f}, sigma = {s}'.format(k=k, f=f, s=sigma)
                    print 'exp(k^2) = {x}'.format(x=np.exp(np.power(k, 2)))

                    raise

        else:

            mu0 = 1 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/3. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/10.

            mu2 = 0.3333333333333333 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/5. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/14.

            mu4 = 0.2 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/7. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/18.

            mu6 = 0.14285714285714285 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/9. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/22.

            mu8 = 0.1111111111111111 - (np.power(f,2)*np.power(k,2)*np.power(sigma,2))/11. + (np.power(f,4)*np.power(k,4)*np.power(sigma,4))/26.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_2(self, k, f, sigma):

        if abs(sigma) > settings.LowSigmaSeries:

            expf = np.exp(-np.power(f,2)*np.power(k,2)*np.power(sigma,2))
            erff = erf(f*k*sigma)

            with np.errstate(over='raise'):

                try:

                    mu0 = expf*(-15/(4.*np.power(f,2)*np.power(k,2)*np.power(sigma,2))) + erff*((-5*np.sqrt(np.pi)*(-3 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(8.*np.power(f,3)*np.power(k,3)*np.power(sigma,3)))

                    mu2 = expf*((-5*(9 + 4*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(8.*np.power(f,4)*np.power(k,4)*np.power(sigma,4))) + erff*((-5*np.sqrt(np.pi)*(-9 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(16.*np.power(f,5)*np.power(k,5)*np.power(sigma,5)))

                    mu4 = expf*((-5*(45 + 24*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 8*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(16.*np.power(f,6)*np.power(k,6)*np.power(sigma,6))) + erff*((-15*np.sqrt(np.pi)*(-15 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(32.*np.power(f,7)*np.power(k,7)*np.power(sigma,7)))

                    mu6 = expf*((-5*(315 + 180*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 64*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 16*np.power(f,6)*np.power(k,6)*np.power(sigma,6)))/(32.*np.power(f,8)*np.power(k,8)*np.power(sigma,8))) + erff*((5*(315*np.sqrt(np.pi) - 30*np.power(f,2)*np.power(k,2)*np.sqrt(np.pi)*np.power(sigma,2)))/(64.*np.power(f,9)*np.power(k,9)*np.power(sigma,9)))

                    mu8 = expf*((-5*(2835 + 1680*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 616*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 160*np.power(f,6)*np.power(k,6)*np.power(sigma,6) + 32*np.power(f,8)*np.power(k,8)*np.power(sigma,8)))/(64.*np.power(f,10)*np.power(k,10)*np.power(sigma,10))) + erff*((-525*np.sqrt(np.pi)*(-27 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(128.*np.power(f,11)*np.power(k,11)*np.power(sigma,11)))

                except FloatingPointError:

                    print 'k = {k}, f = {f}, sigma = {s}'.format(k=k, f=f, s=sigma)
                    print 'exp(k^2) = {x}'.format(x=np.exp(np.power(k, 2)))

                    raise

        else:

            mu0 = (-2*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/3. + (2*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/7.

            mu2 = 0.6666666666666666 - (4*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/7. + (5*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/21.

            mu4 = 0.5714285714285714 - (10*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/21. + (20*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/99.

            mu6 = 0.47619047619047616 - (40*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/99. + (25*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu8 = 0.40404040404040403 - (50*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (2*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/13.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_4(self, k, f, sigma):

        if abs(sigma) > settings.LowSigmaSeries:

            expf = np.exp(-np.power(f,2)*np.power(k,2)*np.power(sigma,2))
            erff = erf(f*k*sigma)

            with np.errstate(over='raise'):

                try:

                    mu0 = expf*((-45*(21 + 2*np.power(f,2)*np.power(k,2)*np.power(sigma,2)))/(32.*np.power(f,4)*np.power(k,4)*np.power(sigma,4))) + erff*((27*np.sqrt(np.pi)*(35 - 20*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(64.*np.power(f,5)*np.power(k,5)*np.power(sigma,5)))

                    mu2 = expf*((-9*(525 + 170*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 32*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(64.*np.power(f,6)*np.power(k,6)*np.power(sigma,6))) + erff*((27*np.sqrt(np.pi)*(175 - 60*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(128.*np.power(f,7)*np.power(k,7)*np.power(sigma,7)))

                    mu4 = expf*((-9*(3675 + 1550*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 416*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 64*np.power(f,6)*np.power(k,6)*np.power(sigma,6)))/(128.*np.power(f,8)*np.power(k,8)*np.power(sigma,8))) + erff*((9*(3675*np.sqrt(np.pi) - 900*np.power(f,2)*np.power(k,2)*np.sqrt(np.pi)*np.power(sigma,2) + 36*np.power(f,4)*np.power(k,4)*np.sqrt(np.pi)*np.power(sigma,4)))/(256.*np.power(f,9)*np.power(k,9)*np.power(sigma,9)))

                    mu6 = expf*((-9*(33075 + 15750*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4800*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 960*np.power(f,6)*np.power(k,6)*np.power(sigma,6) + 128*np.power(f,8)*np.power(k,8)*np.power(sigma,8)))/(256.*np.power(f,10)*np.power(k,10)*np.power(sigma,10))) + erff*((405*np.sqrt(np.pi)*(735 - 140*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(512.*np.power(f,11)*np.power(k,11)*np.power(sigma,11)))

                    mu8 = expf*((-9*(363825 + 185850*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 60480*np.power(f,4)*np.power(k,4)*np.power(sigma,4) + 13440*np.power(f,6)*np.power(k,6)*np.power(sigma,6) + 2176*np.power(f,8)*np.power(k,8)*np.power(sigma,8) + 256*np.power(f,10)*np.power(k,10)*np.power(sigma,10)))/(512.*np.power(f,12)*np.power(k,12)*np.power(sigma,12))) + erff*((2835*np.sqrt(np.pi)*(1155 - 180*np.power(f,2)*np.power(k,2)*np.power(sigma,2) + 4*np.power(f,4)*np.power(k,4)*np.power(sigma,4)))/(1024.*np.power(f,13)*np.power(k,13)*np.power(sigma,13)))

                except FloatingPointError:

                    print 'k = {k}, f = {f}, sigma = {s}'.format(k=k, f=f, s=sigma)
                    print 'exp(k^2) = {x}'.format(x=np.exp(np.power(k, 2)))

                    raise

        else:

            mu0 = (4*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/35.

            mu2 = (-8*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/35. + (12*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/77.

            mu4 = 0.22857142857142856 - (24*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/77. + (24*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu6 = 0.3116883116883117 - (48*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (24*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/143.

            mu8 = 0.3356643356643357 - (48*np.power(f,2)*np.power(k,2)*np.power(sigma,2))/143. + (36*np.power(f,4)*np.power(k,4)*np.power(sigma,4))/221.

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
        
        initial_sigmav = np.array([1.0])
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

