from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft as heft

import WizCOLA
import EFT

import params as lin


class linear_pipeline(heft.cosmosis_pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # set up container of k-samples for WizCOLA
        ks = WizCOLA.ksamples(my_config)

        # build data container
        data = WizCOLA.products(my_config, ks)

        # build theory container
        theory = EFT.EFT_products(my_config, ks)

        # call base class constructor
        super(linear_pipeline, self).__init__(my_name, data, theory)


    def execute(self, block):

        # extract parameters from datablock
        b1 = block['bias_parameters', 'b1']

        params = lin.make_params({'b1': b1})

        return super(linear_pipeline, self).compute(block, params, b1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(linear_pipeline, self).cleanup()


# register this module with the CosmoSIS core
declare_module(linear_pipeline)
