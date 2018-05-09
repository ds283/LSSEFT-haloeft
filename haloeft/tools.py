import numpy as np

import os
import csv

from base import base
from config import settings


class tools(base):

    def __init__(self, model_name, data, theory):

        super(tools, self).__init__(data, theory)

        self.model_name = model_name


    def compute_chisq_variation(self, P0, P2, P4, label_prefix):

        P = np.concatenate((P0, P2, P4))

        # set up empty dictionary to hold return values
        rval = {}

        # set up mask: we will gradually add True values to this as we step through the different
        # k-sample points that contribute to the final chi-square value
        # notice the mask is for just one P_ell individually, not the concatenated group (P0, P2, P4)
        i_mask = np.array([False for i in xrange(len(self.data.k_sample.mean_ks))])

        convolved_Pk = {}
        # cache convolved theory power spectra
        for r in self.data.regions:

            Pconv = np.dot(self.data.convs[r], P)

            if np.any(abs(Pconv) > settings.PconvCeiling):
                print '!! region {tag}: Pconv very large: largest value = {l}, P.max = {P}, conv.max = {c}'.format(
                    tag=r, l=abs(Pconv).max(), P=abs(P).max(), c=abs(self.data.convs[r]).max())
                print Pconv
                raise RuntimeError

            convolved_Pk[r] = Pconv[self.data.k_sample.conv_to_means_mask]

        # step through each available sample point
        for i in xrange(len(i_mask)):

            # is this sample point included in the fit mask?
            i_mask[i] = self.data.k_sample.mean_fit_mask[i]

            # if so, compute the chi-square up to this k-sample point; otherwise, ignore it
            if i_mask[i]:

                # extend mask for the concatenated group (P0, P2, P4)
                i_mask_full = np.concatenate((i_mask, i_mask, i_mask))

                # zero accumulator for chi-square
                chisq = 0.0

                # loop over all data realizations:
                for r in self.data.regions:
                    # extract data products for this region
                    means = self.data.fit_means[r]
                    inv_cov = self.data.fit_inv_covs[r]

                    # cut raw region used in fit (ie. from 0.01 to 0.29) out of full convolved vector
                    Ptheory = convolved_Pk[r]

                    # cut down to just the k-sample points we are using here
                    Ptheory_cut = Ptheory[i_mask_full]

                    # cut a mask for the data means and data covariance matrix out of the mask for the raw region
                    i_mask_means = i_mask_full[self.data.k_sample.mean_fit_mask]
                    means_cut = means[i_mask_means]
                    inv_cov_cut = inv_cov[i_mask_means, :][:, i_mask_means]

                    # take difference between theory and data, and process into a chi-square
                    Delta_cut = Ptheory_cut - means_cut
                    chisq += np.dot(np.dot(Delta_cut, inv_cov_cut), Delta_cut)

                rval[label_prefix + '_' + self.data.k_sample.labels[i]] = chisq

        return rval


    def make_goodness_of_fit_plot(self, plot_file, P0, P2, P4):

        P = np.concatenate((P0, P2, P4))

        # compute total chi-square value
        total_chisq = -2.0 * self.compute_likelihood(P0, P2, P4, type='fit')

        with open(plot_file, 'w') as f:

            print ':: generating theory + data fit comparison file "{f}"'.format(f=plot_file)

            # add labels for global columns
            labels = ['ireal', 'chisq', 'model',
                      'k', 'P0', 'P2', 'P4',
                      'P0_theory_avg', 'P2_theory_avg', 'P4_theory_avg',
                      'P0_WizCOLA_avg', 'P2_WizCOLA_avg', 'P4_WizCOLA_avg',
                      'P0_WizCOLA_err_avg', 'P2_WizCOLA_err_avg', 'P4_WizCOLA_err_avg',
                      'P0_ensemble_avg', 'P2_ensemble_avg', 'P4_ensemble_avg',
                      'P0_ensemble_err_avg', 'P2_ensemble_err_avg', 'P4_ensemble_err_avg',
                      'P0_Delta_avg', 'P2_Delta_avg', 'P4_Delta_avg']

            # supplement with labels for per-region columns
            for r in self.data.regions:
                labels += [r + '_P0_theory', r + '_P2_theory', r + '_P4_theory',
                           r + '_P0_WizCOLA', r + '_P2_WizCOLA', r + '_P4_WizCOLA',
                           r + '_P0_WizCOLA_err', r + '_P2_WizCOLA_err', r + '_P4_WizCOLA_err',
                           r + '_P0_Delta', r + '_P2_Delta', r + '_P4_Delta']

            # open a CSV writer for these columns and emit a header
            writer = csv.DictWriter(f, labels, restval='MISSING')
            writer.writeheader()

            convolved_Pk = {}
            # cache convolved theory power spectra
            for r in self.data.regions:

                Pconv = np.dot(self.data.convs[r], P)

                if np.any(abs(Pconv) > settings.PconvCeiling):
                    print '!! region {tag}: Pconv very large: largest value = {l}, P.max = {P}, conv.max = {c}'.format(
                        tag=r, l=abs(Pconv).max(), P=abs(P).max(), c=abs(self.data.convs[r]).max())
                    print Pconv
                    raise RuntimeError

                convolved_Pk[r] = Pconv

            # loop over each k sample point to be included in the analysis,
            # and then process its contribution for each region
            for i in xrange(len(self.data.k_sample.mean_ks)):

                row = {'k': self.data.k_sample.mean_ks[i], 'P0': P0[i], 'P2': P2[i], 'P4': P4[i]}

                if i == 0:
                    row.update({'ireal': self.data.get_realization(), 'chisq': total_chisq, 'model': self.model_name})
                else:
                    row.update({'ireal': -99, 'chisq': -99, 'model': '-'})

                # initialize accumulators used to sum over all realizations
                P0_theory_tot = 0
                P2_theory_tot = 0
                P4_theory_tot = 0

                P0_WizCOLA_tot = 0
                P2_WizCOLA_tot = 0
                P4_WizCOLA_tot = 0

                P0_ensemble_tot = 0
                P2_ensemble_tot = 0
                P4_ensemble_tot = 0

                P0_WizCOLA_var_tot = 0
                P2_WizCOLA_var_tot = 0
                P4_WizCOLA_var_tot = 0

                P0_ensemble_var_tot = 0
                P2_ensemble_var_tot = 0
                P4_ensemble_var_tot = 0

                P0_Delta_tot = 0
                P2_Delta_tot = 0
                P4_Delta_tot = 0

                # loop over realizations
                for r in self.data.regions:
                    means = self.data.raw_means[r]
                    variances = self.data.raw_variance[r]

                    ensemble_means = self.data.ensemble_means[r]
                    ensemble_variances = self.data.ensemble_variance[r]

                    # read previously cached convolved power spectra
                    Pconv = convolved_Pk[r]

                    # pick out P0, P2, P4 values for this k-sample point
                    P0_theory = Pconv[0 * self.data.k_sample.nbinc + i]
                    P2_theory = Pconv[1 * self.data.k_sample.nbinc + i]
                    P4_theory = Pconv[2 * self.data.k_sample.nbinc + i]

                    P0_data = means[0 * self.data.k_sample.nbin + i]
                    P2_data = means[1 * self.data.k_sample.nbin + i]
                    P4_data = means[2 * self.data.k_sample.nbin + i]

                    P0_ensemble = ensemble_means[0 * self.data.k_sample.nbin + i]
                    P2_ensemble = ensemble_means[1 * self.data.k_sample.nbin + i]
                    P4_ensemble = ensemble_means[2 * self.data.k_sample.nbin + i]

                    P0_data_var = variances[0 * self.data.k_sample.nbin + i]
                    P2_data_var = variances[1 * self.data.k_sample.nbin + i]
                    P4_data_var = variances[2 * self.data.k_sample.nbin + i]

                    P0_ensemble_var = ensemble_variances[0 * self.data.k_sample.nbin + i]
                    P2_ensemble_var = ensemble_variances[1 * self.data.k_sample.nbin + i]
                    P4_ensemble_var = ensemble_variances[2 * self.data.k_sample.nbin + i]

                    row.update({r + '_P0_theory': P0_theory, r + '_P0_WizCOLA': P0_data,
                                r + '_P0_WizCOLA_err': np.sqrt(P0_data_var)})
                    row.update({r + '_P2_theory': P2_theory, r + '_P2_WizCOLA': P2_data,
                                r + '_P2_WizCOLA_err': np.sqrt(P2_data_var)})
                    row.update({r + '_P4_theory': P4_theory, r + '_P4_WizCOLA': P4_data,
                                r + '_P4_WizCOLA_err': np.sqrt(P4_data_var)})

                    P0_delta = P0_data - P0_theory
                    P2_delta = P2_data - P2_theory
                    P4_delta = P4_data - P4_theory

                    row.update({r + '_P0_Delta': P0_delta})
                    row.update({r + '_P2_Delta': P2_delta})
                    row.update({r + '_P4_Delta': P4_delta})

                    P0_theory_tot += P0_theory
                    P2_theory_tot += P2_theory
                    P4_theory_tot += P4_theory

                    P0_WizCOLA_tot += P0_data
                    P2_WizCOLA_tot += P2_data
                    P4_WizCOLA_tot += P4_data

                    P0_ensemble_tot += P0_ensemble
                    P2_ensemble_tot += P2_ensemble
                    P4_ensemble_tot += P4_ensemble

                    P0_WizCOLA_var_tot += P0_data_var
                    P2_WizCOLA_var_tot += P2_data_var
                    P4_WizCOLA_var_tot += P4_data_var

                    P0_ensemble_var_tot += P0_ensemble_var
                    P2_ensemble_var_tot += P2_ensemble_var
                    P4_ensemble_var_tot += P4_ensemble_var

                    P0_Delta_tot += P0_delta
                    P2_Delta_tot += P2_delta
                    P4_Delta_tot += P4_delta

                regions = len(self.data.regions)

                row.update({'P0_theory_avg': P0_theory_tot / regions,
                            'P2_theory_avg': P2_theory_tot / regions,
                            'P4_theory_avg': P4_theory_tot / regions,
                            'P0_WizCOLA_avg': P0_WizCOLA_tot / regions,
                            'P2_WizCOLA_avg': P2_WizCOLA_tot / regions,
                            'P4_WizCOLA_avg': P4_WizCOLA_tot / regions,
                            'P0_WizCOLA_err_avg': np.sqrt(P0_WizCOLA_var_tot) / regions,
                            'P2_WizCOLA_err_avg': np.sqrt(P2_WizCOLA_var_tot) / regions,
                            'P4_WizCOLA_err_avg': np.sqrt(P4_WizCOLA_var_tot) / regions,
                            'P0_ensemble_avg': P0_ensemble_tot / regions,
                            'P2_ensemble_avg': P2_ensemble_tot / regions,
                            'P4_ensemble_avg': P4_ensemble_tot / regions,
                            'P0_ensemble_err_avg': np.sqrt(P0_ensemble_var_tot) / regions,
                            'P2_ensemble_err_avg': np.sqrt(P2_ensemble_var_tot) / regions,
                            'P4_ensemble_err_avg': np.sqrt(P4_ensemble_var_tot) / regions,
                            'P0_Delta_avg': P0_Delta_tot / regions,
                            'P2_Delta_avg': P2_Delta_tot / regions,
                            'P4_Delta_avg': P4_Delta_tot / regions})

                writer.writerow(row)
