from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft as heft

import WizCOLA
import EFT

import params as fp


class full_pipeline(heft.cosmosis_pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # set up container of k-samples for WizCOLA
        ks = WizCOLA.ksamples(my_config)

        # build data container
        data = WizCOLA.products(my_config, ks)

        # build theory container
        theory = EFT.EFT_products(my_config, ks)

        # call base class constructor
        super(full_pipeline, self).__init__(my_name, data, theory)


    def execute(self, block):

        # extract parameters from datablock
        b1_1 = block['bias_parameters', 'b1_1']
        b1_2 = block['bias_parameters', 'b1_2']
        b1_3 = block['bias_parameters', 'b1_3']
        b2_2 = block['bias_parameters', 'b2_2']
        bG2_2 = block['bias_parameters', 'bG2_2']
        bG2_3 = block['bias_parameters', 'bG2_3']

        # build coefficient dictionary
        params = fp.make_params({'b1_1': b1_1, 'b1_2': b1_2, 'b1_3': b1_3, 'b2_2': b2_2, 'bG2_2': bG2_2, 'bG2_3': bG2_3})

        return super(full_pipeline, self).compute(block, params, b1_1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(full_pipeline, self).cleanup()


# register this module with the CosmoSIS core
declare_module(full_pipeline)
