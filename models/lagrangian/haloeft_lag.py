from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft as heft

import WizCOLA
import EFT

import lag_params as lag


class HaloEFT(heft.cosmosis_pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # set up container of k-samples for WizCOLA
        ks = WizCOLA.ksamples(my_config)

        # build data container
        data = WizCOLA.products(my_config, ks)

        # build theory container
        theory = EFT.EFT_products(my_config, ks)

        # call base class constructor
        super(HaloEFT, self).__init__(my_name, data, theory)


    def execute(self, block):

        # extract parameters from datablock
        b1 = block['bias_parameters', 'b1']
        b2 = block['bias_parameters', 'b2']
        bs2 = block['bias_parameters', 'bs2']
        b3nl = block['bias_parameters', 'b3nl']

        params = lag.make_params({'b1': b1, 'b2': b2, 'bs2': bs2, 'b3nl': b3nl})

        return super(HaloEFT, self).compute(block, params, b1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(HaloEFT, self).cleanup()


# register this module with the CosmoSIS core
declare_module(HaloEFT)
