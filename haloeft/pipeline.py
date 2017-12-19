from base import base

class cosmosis_pipeline(base):

    def __init__(self, my_name, data, theory):

        super(cosmosis_pipeline, self).__init__(data, theory)

        self.__mod_name = my_name


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


    def cleanup(self):

        return 0
