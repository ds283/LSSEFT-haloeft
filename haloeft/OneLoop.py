import pipeline as cp

import WizCOLA
import RedshiftModels.OneLoop as ZhengSong


class pipeline(cp.cosmosis_pipeline):

    def __init__(self, my_config, my_name):

        # set up container of k-samples for WizCOLA
        ks = WizCOLA.ksamples(my_config)

        # build data container
        data = WizCOLA.products(my_config, ks)

        # build theory container
        theory = OneLoop.theory(my_config, ks)

        # call base class constructor
        super(pipeline, self).__init__(my_name, data, theory)


    def cleanup(self):

        # pass on to base class cleanup
        super(pipeline, self).cleanup()
