import numpy as np


# container class for k-sample points
class ksamples(object):

    def __init__(self, my_config):

        self.nbin = 15
        self.nbinc = 25

        self.WiggleZ_mean_ks = np.linspace(0.01, 0.29, self.nbin)
        self.WiggleZ_conv_ks = np.linspace(0.01, 0.49, self.nbinc)

        self.labels = ['0.01', '0.03', '0.05', '0.07', '0.09',
                       '0.11', '0.13', '0.15', '0.17', '0.19',
                       '0.21', '0.23', '0.25', '0.27', '0.29']


        # BUILD AND CACHE MASKS TO SELECT DIFFERENT PARTS OF THE DATA PRODUCTS

        fit_kmin = my_config["HaloEFT", "fit_kmin"]
        fit_kmax = my_config["HaloEFT", "fit_kmax"]

        mam = np.all([self.WiggleZ_mean_ks > fit_kmin, self.WiggleZ_mean_ks <= fit_kmax], axis=0)
        cam = np.all([self.WiggleZ_conv_ks > fit_kmin, self.WiggleZ_conv_ks <= fit_kmax], axis=0)
        self.mean_fit_mask = np.concatenate((mam, mam, mam))
        self.conv_fit_mask = np.concatenate((cam, cam, cam))

        ren_kmin = my_config["HaloEFT", "renormalize_kmin"]
        ren_kmax = my_config["HaloEFT", "renormalize_kmax"]

        mrm = np.all([self.WiggleZ_mean_ks > ren_kmin, self.WiggleZ_mean_ks <= ren_kmax], axis=0)
        crm = np.all([self.WiggleZ_conv_ks > ren_kmin, self.WiggleZ_conv_ks <= ren_kmax], axis=0)
        self.mean_ren_mask = np.concatenate((mrm, mrm, mrm))
        self.conv_ren_mask = np.concatenate((crm, crm, crm))

        cmm = self.WiggleZ_conv_ks <= 0.30
        self.conv_to_means_mask = np.concatenate((cmm, cmm, cmm))
