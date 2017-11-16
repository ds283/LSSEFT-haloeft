from cosmosis.runtime.declare import declare_module
import haloeft as heft
import full_params as fp


class HaloEFT(heft.HaloEFT_core):

    def __init__(self, my_config, my_name):

        # call base class constructor
        super(HaloEFT, self).__init__(my_config, my_name)


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

        return super(HaloEFT, self).compute(block, params, b1_1)


    def cleanup(self):

        # pass on to base class cleanup
        super(HaloEFT, self).cleanup()


# register this module with the CosmoSIS core
declare_module(HaloEFT)
