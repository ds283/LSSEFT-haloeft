from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft.Kaiser as Kaiser

import models.NonlocalAll.params as fp


class pipeline(Kaiser.pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # call base class constructor
        super(pipeline, self).__init__(my_config, my_name)


    def execute(self, block):

        # extract parameters from datablock
        b1_1 = block['bias_parameters', 'b1_1']
        b1_2 = block['bias_parameters', 'b1_2']
        b1_3 = block['bias_parameters', 'b1_3']
        b2_2 = block['bias_parameters', 'b2_2']
        b2_3 = block['bias_parameters', 'b2_3']
        b3 = block['bias_parameters', 'b3']
        bG2_2 = block['bias_parameters', 'bG2_2']
        bG2_3 = block['bias_parameters', 'bG2_3']
        bdG2 = block['bias_parameters', 'bdG2']
        bGamma3 = block['bias_parameters', 'bGamma3']

        # build coefficient dictionary
        params = fp.make_params({'b1_1': b1_1, 'b1_2': b1_2, 'b1_3': b1_3, 'b2_2': b2_2, 'b2_3': b2_3,
                                 'b3': b3, 'bG2_2': bG2_2, 'bG2_3': bG2_3, 'bdG2': bdG2, 'bGamma3': bGamma3})

        return super(pipeline, self).compute(block, params, b1_1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(pipeline, self).cleanup()


# register this module with the CosmoSIS core
declare_module(pipeline)
