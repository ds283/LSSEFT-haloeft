import numpy as np
from astropy.io import ascii

from config import settings
from utils import matrix as mtx


# container class for WizCOLA data products
class products(object):

    def __init__(self, my_config, k_sample):

        self.k_sample = k_sample

        # extract desired realization from datablock
        self.__realization = my_config["HaloEFT", "realization"]


        # READ AND CACHE WIGGLEZ DATA PRODUCTS

        h_means = {}
        h_cov = {}

        h_means.update({'1h': my_config["HaloEFT", "h1_means"]})
        h_cov.update({'1h': my_config["HaloEFT", "h1_matrix"]})

        h_means.update({'3h': my_config["HaloEFT", "h3_means"]})
        h_cov.update({'3h': my_config["HaloEFT", "h3_matrix"]})

        h_means.update({'9h': my_config["HaloEFT", "h9_means"]})
        h_cov.update({'9h': my_config["HaloEFT", "h9_matrix"]})

        h_means.update({'11h': my_config["HaloEFT", "h11_means"]})
        h_cov.update({'11h': my_config["HaloEFT", "h11_matrix"]})

        h_means.update({'15h': my_config["HaloEFT", "h15_means"]})
        h_cov.update({'15h': my_config["HaloEFT", "h15_matrix"]})

        h_means.update({'22h': my_config["HaloEFT", "h22_means"]})
        h_cov.update({'22h': my_config["HaloEFT", "h22_matrix"]})

        # object-local storage for data products, held as dictionaries indexed by region label

        # 1. Power spectrum mean values: (a) raw, (b) cut down to region needed for likelihood fit,
        # (c) cut down to region needed for renormalization
        self.raw_means = {}
        self.fit_means = {}
        self.ren_means = {}

        # 2. Inverse covariance matrices: (a) cut down to region needed for likelihood fit,
        # (b) cut down to region needed for renormalization
        self.fit_inv_covs = {}
        self.ren_inv_covs = {}

        # 4. Measurement errors on the power spectrum mean value
        self.raw_variance = {}

        # 5. Ensemble average measurements, per region
        self.ensemble_means = {}
        self.ensemble_variance = {}

        # 6. Convolution matrix taking theory outputs to each WizCOLA region
        # (eg. accounts for geometry and selection function)
        self.convs = {}

        # 7. Region labels
        self.regions = ['1h', '3h', '9h', '11h', '15h', '22h']

        # perform import
        for r in self.regions:

            self.__import_WiggleZ_data(r, h_means[r], h_cov[r], my_config)


    def __import_WiggleZ_data(self, tag, means_file, matrix_file, my_config):

        # number of bins in P_ell measurements and covariance matrix
        nbin = self.k_sample.nbin

        # number of bins in convolution  matrix
        nbinc = self.k_sample.nbinc

        # read in the WizCOLA data for this region as astropy.table.Table instances
        mean_table = ascii.read(means_file, Reader=ascii.NoHeader,
                                names=['ireal', 'k', 'P0', 'P0err', 'P2', 'P2err', 'P4', 'P4err'])

        cov_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=0, data_end=(3*nbin)**2,
                               names=['i', 'j', 'value'])

        conv_table = ascii.read(matrix_file, Reader=ascii.NoHeader, data_start=(3*nbin)**2, data_end=(3*nbin)**2+(3*nbinc)**2,
                                names=['i', 'j', 'value'])

        # want to extract the power spectrum for a given realization, so group mean_table by realization
        mean_by_ireal = mean_table.group_by('ireal')

        # get list of all realizations
        all_realizations = mean_by_ireal.groups.keys['ireal']

        # use mask to extract Pk means for just this realization
        mean_this_realization = mean_by_ireal.groups[self.__realization]

        # sort into ascending order of k
        mean_this_realization.sort('k')

        # extract P0, P2, P4 columns as numpy arrays
        this_P0 = np.asarray(mean_this_realization['P0'])
        this_P2 = np.asarray(mean_this_realization['P2'])
        this_P4 = np.asarray(mean_this_realization['P4'])

        # concatenate these arrays to form a single block consisting of P0 values followed by P2, P4 values
        Pblock = np.concatenate((this_P0, this_P2, this_P4))

        # extract relevant components -- 'raw' is everything, but 'fit' and 'ren' use only a subset of values
        self.raw_means[tag] = Pblock
        self.fit_means[tag] = Pblock[self.k_sample.mean_fit_mask]
        self.ren_means[tag] = Pblock[self.k_sample.mean_ren_mask]


        # construct empty matrices to hold covariance and convolution matrices after read-in
        covariance = np.empty((3*nbin, 3*nbin))
        convolution = np.empty((3*nbinc, 3*nbinc))

        # import matrices
        mtx.import_matrix(covariance, cov_table, tag, 'covariance')
        mtx.import_matrix(convolution, conv_table, tag, 'convolution')


        # strip out diagonal entries of covariance matrix as error estimates
        variances = np.diagonal(covariance)
        self.raw_variance[tag] = variances


        # cut covariance matrix down only to those values used in the likelihood fit
        # (notice we do this *before* inversion)
        covariance_cut_fit = (covariance[self.k_sample.mean_fit_mask, :])[:, self.k_sample.mean_fit_mask]

        # invert and store
        self.fit_inv_covs[tag] = mtx.invert_covariance_matrix(covariance_cut_fit, tag, 'all')


        # cut covariance matrix down only to those values used in renormalization
        # (notice we do this *before* inversion)
        covariance_cut_ren = (covariance[self.k_sample.mean_ren_mask, :])[:, self.k_sample.mean_ren_mask]

        # invert and store
        self.ren_inv_covs[tag] = mtx.invert_covariance_matrix(covariance_cut_ren, tag, 'ren')

        # store convolution matrix, first checking whether it contains any anomalously-large elements
        # that would be symptomatic of an error
        if np.any(abs(convolution) > settings.ConvCeiling):

            print '!! region {tag}: convolution matrix contains very large element = {v}'.format(tag=tag, v=abs(convolution).max())
            raise RuntimeError

        self.convs[tag] = convolution


        # extract ensemble average power spectrum for this region, which is stored as realization 0
        ensemble_means = mean_by_ireal.groups[0]
        ensemble_means.sort('k')

        ensemble_P0 = np.asarray(ensemble_means['P0'])
        ensemble_P2 = np.asarray(ensemble_means['P2'])
        ensemble_P4 = np.asarray(ensemble_means['P4'])

        ensemble_Pk = np.concatenate((ensemble_P0, ensemble_P2, ensemble_P4))

        self.ensemble_means[tag] = ensemble_Pk
        self.ensemble_variance[tag] = variances / (len(all_realizations) - 1.0)


    def get_realization(self):

        return self.__realization

