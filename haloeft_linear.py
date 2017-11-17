from cosmosis.runtime.declare import declare_module
import haloeft as heft
import linear_params as lin


class HaloEFT(heft.HaloEFT_core):

    def __init__(self, my_config, my_name):

        # call base class constructor
        super(HaloEFT, self).__init__(my_config, my_name)


    def execute(self, block):

        # extract parameters from datablock
        b1 = block['bias_parameters', 'b1']

        params = lin.make_params({'b1': b1})

        return super(HaloEFT, self).compute(block, params, b1)


    def cleanup(self):

        # pass on to base class cleanup
        super(HaloEFT, self).cleanup()


# register this module with the CosmoSIS core
declare_module(HaloEFT)
