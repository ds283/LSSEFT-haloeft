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


    def __ell_0(self, x):

        if abs(x) > settings.RSDSeriesSwitchover:

            expx = np.exp(-x * x)
            erfx = erf(x)

            with np.errstate(over='raise'):

                try:

                    mu0 = expx * (0) + erfx * (np.sqrt(np.pi)/(2.*x))

                    mu2 = expx * (-1/(2.*np.power(x,2))) + erfx * (np.sqrt(np.pi)/(4.*np.power(x,3)))

                    mu4 = expx * (-(3 + 2*np.power(x,2))/(4.*np.power(x,4))) + erfx * ((3*np.sqrt(np.pi))/(8.*np.power(x,5)))

                    mu6 = expx * (-(15 + 10*np.power(x,2) + 4*np.power(x,4))/(8.*np.power(x,6))) + erfx * ((15*np.sqrt(np.pi))/(16.*np.power(x,7)))

                    mu8 = expx * (-(105 + 70*np.power(x,2) + 28*np.power(x,4) + 8*np.power(x,6))/(16.*np.power(x,8))) + erfx * ((105*np.sqrt(np.pi))/(32.*np.power(x,9)))

                except FloatingPointError:

                    print 'x = {x}'.format(x=x)
                    print 'exp(-x^2) = {x}'.format(x=expx)
                    print 'erf(x) = {x}'.format(x=erfx)

                    raise

        else:

            mu0 = 1 - np.power(x,2)/3. + np.power(x,4)/10. - np.power(x,6)/42. + np.power(x,8)/216. - np.power(x,10)/1320.

            mu2 = 0.3333333333333333 - np.power(x,2)/5. + np.power(x,4)/14. - np.power(x,6)/54. + np.power(x,8)/264. - np.power(x,10)/1560.

            mu4 = 0.2 - np.power(x,2)/7. + np.power(x,4)/18. - np.power(x,6)/66. + np.power(x,8)/312. - np.power(x,10)/1800.

            mu6 = 0.14285714285714285 - np.power(x,2)/9. + np.power(x,4)/22. - np.power(x,6)/78. + np.power(x,8)/360. - np.power(x,10)/2040.

            mu8 = 0.1111111111111111 - np.power(x,2)/11. + np.power(x,4)/26. - np.power(x,6)/90. + np.power(x,8)/408. - np.power(x,10)/2280.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_2(self, x):

        if abs(x) > settings.RSDSeriesSwitchover:

            expx = np.exp(-x * x)
            erfx = erf(x)

            with np.errstate(over='raise'):

                try:

                    mu0 = expx * (-15/(4.*np.power(x,2))) + erfx * ((-5*np.sqrt(np.pi)*(-3 + 2*np.power(x,2)))/(8.*np.power(x,3)))

                    mu2 = expx * ((-5*(9 + 4*np.power(x,2)))/(8.*np.power(x,4))) + erfx * ((-5*np.sqrt(np.pi)*(-9 + 2*np.power(x,2)))/(16.*np.power(x,5)))

                    mu4 = expx * ((-5*(45 + 24*np.power(x,2) + 8*np.power(x,4)))/(16.*np.power(x,6))) + erfx * ((-15*np.sqrt(np.pi)*(-15 + 2*np.power(x,2)))/(32.*np.power(x,7)))

                    mu6 = expx * ((-5*(315 + 180*np.power(x,2) + 64*np.power(x,4) + 16*np.power(x,6)))/(32.*np.power(x,8))) + erfx * ((5*(315*np.sqrt(np.pi) - 30*np.sqrt(np.pi)*np.power(x,2)))/(64.*np.power(x,9)))

                    mu8 = expx * ((-5*(2835 + 1680*np.power(x,2) + 616*np.power(x,4) + 160*np.power(x,6) + 32*np.power(x,8)))/(64.*np.power(x,10))) + erfx * ((5*(2835*np.sqrt(np.pi) - 210*np.sqrt(np.pi)*np.power(x,2)))/(128.*np.power(x,11)))

                except FloatingPointError:

                    print 'x = {x}'.format(x=x)
                    print 'exp(-x^2) = {x}'.format(x=expx)
                    print 'erf(x) = {x}'.format(x=erfx)

                    raise

        else:

            mu0 = (-2*np.power(x,2))/3. + (2*np.power(x,4))/7. - (5*np.power(x,6))/63. + (5*np.power(x,8))/297. - (5*np.power(x,10))/1716.

            mu2 = 0.6666666666666666 - (4*np.power(x,2))/7. + (5*np.power(x,4))/21. - (20*np.power(x,6))/297. + (25*np.power(x,8))/1716. - np.power(x,10)/390.

            mu4 = 0.5714285714285714 - (10*np.power(x,2))/21. + (20*np.power(x,4))/99. - (25*np.power(x,6))/429. + np.power(x,8)/78. - (7*np.power(x,10))/3060.

            mu6 = 0.47619047619047616 - (40*np.power(x,2))/99. + (25*np.power(x,4))/143. - (2*np.power(x,6))/39. + (7*np.power(x,8))/612. - (2*np.power(x,10))/969.

            mu8 = 0.40404040404040403 - (50*np.power(x,2))/143. + (2*np.power(x,4))/13. - (7*np.power(x,6))/153. + (10*np.power(x,8))/969. - np.power(x,10)/532.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __ell_4(self, x):

        if abs(x) > settings.RSDSeriesSwitchover:

            expx = np.exp(-x * x)
            erfx = erf(x)

            with np.errstate(over='raise'):

                try:

                    mu0 = expx * ((-45*(21 + 2*np.power(x,2)))/(32.*np.power(x,4))) + erfx * ((27*np.sqrt(np.pi)*(35 - 20*np.power(x,2) + 4*np.power(x,4)))/(64.*np.power(x,5)))

                    mu2 = expx * ((-9*(525 + 170*np.power(x,2) + 32*np.power(x,4)))/(64.*np.power(x,6))) + erfx * ((27*np.sqrt(np.pi)*(175 - 60*np.power(x,2) + 4*np.power(x,4)))/(128.*np.power(x,7)))

                    mu4 = expx * ((-9*(3675 + 1550*np.power(x,2) + 416*np.power(x,4) + 64*np.power(x,6)))/(128.*np.power(x,8))) + erfx * ((9*(3675*np.sqrt(np.pi) - 900*np.sqrt(np.pi)*np.power(x,2) + 36*np.sqrt(np.pi)*np.power(x,4)))/(256.*np.power(x,9)))

                    mu6 = expx * ((-9*(33075 + 15750*np.power(x,2) + 4800*np.power(x,4) + 960*np.power(x,6) + 128*np.power(x,8)))/(256.*np.power(x,10))) + erfx * ((405*np.sqrt(np.pi)*(735 - 140*np.power(x,2) + 4*np.power(x,4)))/(512.*np.power(x,11)))

                    mu8 = expx * ((-9*(363825 + 185850*np.power(x,2) + 60480*np.power(x,4) + 13440*np.power(x,6) + 2176*np.power(x,8) + 256*np.power(x,10)))/(512.*np.power(x,12))) + erfx * ((2835*np.sqrt(np.pi)*(1155 - 180*np.power(x,2) + 4*np.power(x,4)))/(1024.*np.power(x,13)))

                except FloatingPointError:

                    print 'x = {x}'.format(x=x)
                    print 'exp(-x^2) = {x}'.format(x=expx)
                    print 'erf(x) = {x}'.format(x=erfx)

                    raise

        else:

            mu0 = (4*np.power(x,4))/35. - (4*np.power(x,6))/77. + (2*np.power(x,8))/143. - (2*np.power(x,10))/715.

            mu2 = (-8*np.power(x,2))/35. + (12*np.power(x,4))/77. - (8*np.power(x,6))/143. + (2*np.power(x,8))/143. - (3*np.power(x,10))/1105.

            mu4 = 0.22857142857142856 - (24*np.power(x,2))/77. + (24*np.power(x,4))/143. - (8*np.power(x,6))/143. + (3*np.power(x,8))/221. - (21*np.power(x,10))/8075.

            mu6 = 0.3116883116883117 - (48*np.power(x,2))/143. + (24*np.power(x,4))/143. - (12*np.power(x,6))/221. + (21*np.power(x,8))/1615. - (4*np.power(x,10))/1615.

            mu8 = 0.3356643356643357 - (48*np.power(x,2))/143. + (36*np.power(x,4))/221. - (84*np.power(x,6))/1615. + (4*np.power(x,8))/323. - (36*np.power(x,10))/15295.

        return np.array([mu0, mu2, mu4, mu6, mu8])


    def __dotp(self, a, b):

        # a and b are expected to be matrices, so the produce a[i]*b[i] will extract the ith row from each
        # matrix and perform element-by-element multiplication along the row.
        # Then we sum over rows.

        zip = [ a[i]*b[i] for i in xrange(len(a))]

        P = zip[0].copy()
        for a in zip[1:]:
            P += a

        return P


    def build_theory_P_ell(self, coeffs, values, blinear):

        f = self.data.f
        sigmav = values['sigmav']
        
        ks = self.data.k_sample.conv_ks

        # cache expensive transcendental functions
        xs = ks * f * sigmav

        # First we assemble a matrix of transformation coefficients f_{n,ell} for each k sample point,
        # translated to x = k f sigma
        # The coefficients are defined by
        #    f_{n, ell}(x) = (2ell + 1)/2 Integrate[exp(-x^2 mu^2) mu^n LegP(ell, mu)]

        ell0_coeff = np.transpose(np.array([self.__ell_0(x) for x in xs]))
        ell2_coeff = np.transpose(np.array([self.__ell_2(x) for x in xs]))
        ell4_coeff = np.transpose(np.array([self.__ell_4(x) for x in xs]))

        # Each of these (_after_ the transposition operation) is a numpy array of the form, eg. for ell0_coeff
        # f_{0,0}(x_1)  f_{0,0}(x_2)  f_{0,0}(x_3)  ...  f_{0,0}(x_n)
        # f_{2,0}(x_1)  f_{2,0}(x_2)  f_{2,0}(x_3)  ...  f_{2,0}(x_n)
        # ...
        # f_{8,0}(x_1)  f_{8,0}(x_2)  f_{8,0}(x_3)  ...  f_{8,0}(x_n)
        # where x_1, x_2, ..., x_n are the k-sample points expressed in terms of x.

        # next we construct lists of terms to be added to produce the final multipole power spectra
        # each term coeffs[key] * data is a numpy array giving the powers of mu for a specific
        # bias label at all k sample points, of the form
        # mu0(x_1)  mu0(x_2)  mu0(x_3)  ...  mu0(x_n)
        # mu2(x_1)  mu2(x_2)  mu2(x_3)  ...  mu2(x_n)
        # ...
        # mu8(x_1)  mu8(x_2)  mu8(x_3)  ...  mu8(x_n)

        zip0 = [self.__dotp(ell0_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]
        zip2 = [self.__dotp(ell2_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]
        zip4 = [self.__dotp(ell4_coeff, coeffs[key] * data) for key, data in self.data.payload.iteritems()]

        # the __dotp() function performs a dot product along the first axis of each argument, which is the
        # row axis of each matrix
        # So the dot product will produce a sum of terms of the form
        #   [ f_{0,0}(x_1) mu0(x_1)  f_{0,0}(x_2) mu0(x_2)  f_{0,0}(x_3) mu0(x_3)  ...  f_{0,0}(x_n) mu0(x_n) ]
        #   +
        #   [ f_{2,0}(x_1) mu2(x_1)  f_{2,0}(x_2) mu2(x_2)  f_{2,0}(x_3) mu2(x_3)  ...  f_{2,0}(x_n) mu2(x_n) ]
        #   + ... +
        #   [ f_{8,0}(x_1) mu8(x_1)  f_{8,0}(x_2) mu8(x_2)  f_{8,0}(x_3) mu8(x_3)  ...  f_{8,0}(x_n) mu8(x_n) ]
        # after summing these we get a single vector containing the elements
        #   [ P_0(x_1)  P_0(x_2)  P_0(x_3)  ...  P_0(x_n) ]
        # and likewise for the other multipoles

        # Finally, we sum over all the terms produced in this way, which is a sum over all the different
        # bias coefficients that can appear in the power spectra

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

