from cosmosis.runtime.declare import declare_module
from cosmosis.datablock import names as section_names

import haloeft.ZhengSong as ZhengSong

import models.linear.params as lin


class pipeline(ZhengSong.pipeline):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        # call base class constructor
        super(pipeline, self).__init__(my_config, my_name)


    def execute(self, block):

        # extract parameters from datablock
        b1 = block['bias_parameters', 'b1']

        params = lin.make_params({'b1': b1})

        return super(pipeline, self).compute(block, params, b1, self.likes)


    def cleanup(self):

        # pass on to base class cleanup
        super(pipeline, self).cleanup()


# register this module with the CosmoSIS core
declare_module(pipeline)
