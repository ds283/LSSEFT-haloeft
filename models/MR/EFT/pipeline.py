from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft.EFT as EFT

import models.MR.params as MR


class pipeline(EFT.pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # call base class constructor
        super(pipeline, self).__init__(my_config, my_name)


    def execute(self, block):

        # extract parameters from datablock
        b1 = block['bias_parameters', 'b1']
        b2 = block['bias_parameters', 'b2']
        bs2 = block['bias_parameters', 'bs2']
        b3nl = block['bias_parameters', 'b3nl']

        params = MR.make_params({'b1': b1, 'b2': b2, 'bs2': bs2, 'b3nl': b3nl})

        return super(pipeline, self).compute(block, params, b1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(pipeline, self).cleanup()


# register this module with the CosmoSIS core
declare_module(pipeline)
