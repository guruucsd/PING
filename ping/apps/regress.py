"""
Run the data browser tools from the command-line.
"""
import datetime
import json
import md5
import os
import re

import numpy as np
import requests
import StringIO
from matplotlib import pyplot as plt
from scipy.stats import linregress

from ..access import PINGSession
from ..asymmetry import compute_all_asymmetries


class PINGDataSession(PINGSession):

    def __init__(self, *args, **kwargs):
        super(PINGDataSession, self).__init__(*args, **kwargs)
        self.result_ids = None  # current dictionary of result IDs

    def expert_mode_script(self, X, Y, covariates=[], limits=[]):
        """Produce the necessary 'expert' script"""
        limits = '\n'.join(['data <- data[data$%s,]' % lim for lim in limits])
        return """dependent.measure     = "{Y}"
covariates.usr        = "{covariates}"
covariates.ses        = "FDH_Highest_Education + FDH_3_Household_Income"
covariates.dev        = "DeviceSerialNumber"
covariates.gaf        = "GAF_africa + GAF_amerind + GAF_eastAsia + GAF_oceania + GAF_centralAsia"
independent.variable  = "{X}"
smoothing.interaction = ""

{limits}""".format(X=X, Y=Y, covariates='+'.join(covariates), limits=limits)

    def regress(self, X, Y, covariates=[], limits=[], plot=False, cache_dir='.', force=False, abort_if_done=True):
        """Do the regression remotely (via R), compile result, plot and save locally.
        """
        cookie = np.random.randint(1000)
        payload = {
            '_v': '',
            'cookie': cookie,
            'user_name': self.username,
            'project_name': 'PING',
            'command': '+'.join(covariates),
            'yvalue': Y,
            'functionOf': X,
            'interaction': '',
            'expert': self.expert_mode_script(X, Y, covariates, limits)}
        out_files = [os.path.join(cache_dir, '%s_%s.txt') % (
                         md5.md5(payload['expert']).hexdigest(), ftype)
                     for ftype in ['executeR', 'tsv']]

        # First, generate the regression and store the result.
        if os.path.exists(out_files[0]) and not force:
            if abort_if_done:
                return
            with open(out_files[0], 'rb') as fp:
                r_text = '\n'.join(fp.readlines())
        else:
            url = 'https://ping-dataportal.ucsd.edu/applications/DataExploration/executeR.php'
            self.log("Computing regression for %s vs. %s..." % (X, Y))
            resp = self.sess.post(url, data=payload)
            r_text = str(resp.text)
            with open(out_files[0], 'wb') as fp:
                fp.write(r_text)

        # Parse the regression result
        error_prog = re.compile('.*Error in (.+)<br>Execution halted<br>',
                                re.MULTILINE | re.IGNORECASE)
        good_prog = re.compile('.*<br>full\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)<br>',
                               re.MULTILINE | re.IGNORECASE)
        error_result = error_prog.match(r_text)
        good_result = good_prog.match(r_text)

        if error_result is not None:
            raise Exception(error_result.groups()[0])
        elif good_result is not None:
            outputs = good_result.groups()  # order: n, _, dev, rsq, aic, _, pval
        else:
            raise Exception("Could not parse output: %s" % r_text)

        # Next, retrieve the raw data.
        if os.path.exists(out_files[1]) and not force:
           with open(out_files[1], 'rb') as fp:
                tsv_text = '\n'.join(fp.readlines())
        else:
            self.log('Retrieving raw data...')
            resp = self.sess.get('https://ping-dataportal.ucsd.edu/applications/DataExploration/curves/%s_PING_curves/%s_PING_Corrected%d.tsv?cache=false' % (
                self.username, self.username, cookie))
            tsv_text = str(resp.text)
            with open(out_files[1], 'wb') as fp:
                fp.write(tsv_text)

        # Parse the raw data
        try:
            rec = np.recfromcsv(StringIO.StringIO(tsv_text))
        except Exception as e:
            raise Exception("Failed to parse data response (%s): %s" % (tsv_text, e))
        else:
            self.log("Retrieved: %s" % str(rec.dtype.names))

        # Display the result
        if plot:
            if isinstance(plot, bool):
                ax = plt.figure().gca()
            else:
                ax = plot
            find_one_relationship(rec, X.lower().replace('.', ''), Y.lower().replace('.', ''), plot=ax)
            ax.set_title('%s\n%s' % (ax.get_title(), '\n'.join(covariates[3:])))
        return outputs, rec

    def AI2flds(self, AI):
        if 'MRI_cort_area_ctx' in AI:
            area_name = AI[18:-3]
            return ['MRI_cort_area.ctx.%s.%s' % (hemi, area_name)
                    for hemi in ['lh', 'rh']]
        elif 'MRI_cort_thick_ctx' in AI:
            area_name = AI[19:-3]
            return ['MRI_cort_thick.ctx.%s.%s' % (hemi, area_name)
                    for hemi in ['lh', 'rh']]
        elif 'MRI_subcort_vol' in AI:
            area_name = AI[16:-3].replace('_', '.')
            return ['MRI_subcort_vol.%s.%s' % (hemi, area_name)
                    for hemi in ['Left', 'Right']]
        elif 'DTI_fiber_vol' in AI:
            area_name = AI[14:-3]
            return ['DTI_fiber_vol.%s_%s' % (hemi, area_name)
                    for hemi in ['L', 'R']]
        else:
            raise NotImplementedError(AI)

    def regress_multistep(self, X, Y, covariates=['Age_At_IMGExam', 'Gender', 'FDH_23_Handedness_Prtcpnt'], limits=[], plot=False, cache_dir='download'):
        self.log('Multi-step regression for %s vs. %s' % (X, Y))
        out = []
        # Get relevant covariates
        added_covariates_X = self.AI2flds(X)
        added_covariates_Y = self.AI2flds(Y)

        def add_subplot(fh, *args):
            return fh.add_subplot(*args) if plot else None

        # First, regress with covariates
        fh1 = plt.figure(figsize=(18, 10)) if plot else None
        out.append(self.regress(X, Y, covariates, limits, plot=add_subplot(fh1, 2, 4, 1), cache_dir=cache_dir))
        out.append(self.regress(Y, X, covariates, limits, plot=add_subplot(fh1, 2, 4, 5), cache_dir=cache_dir))
        out.append(self.regress(X, Y, covariates + added_covariates_X, limits, plot=add_subplot(fh1, 2, 4, 2), cache_dir=cache_dir))
        out.append(self.regress(Y, X, covariates + added_covariates_X, limits, plot=add_subplot(fh1, 2, 4, 6), cache_dir=cache_dir))
        out.append(self.regress(X, Y, covariates + added_covariates_Y, limits, plot=add_subplot(fh1, 2, 4, 3), cache_dir=cache_dir))
        out.append(self.regress(Y, X, covariates + added_covariates_Y, limits, plot=add_subplot(fh1, 2, 4, 7), cache_dir=cache_dir))
        out.append(self.regress(X, Y, covariates + added_covariates_X + added_covariates_Y , limits, plot=add_subplot(fh1, 2, 4, 4), cache_dir=cache_dir))
        out.append(self.regress(Y, X, covariates + added_covariates_X + added_covariates_Y , limits, plot=add_subplot(fh1, 2, 4, 8), cache_dir=cache_dir))

        # Then regress on each side
        fh2 = plt.figure(figsize=(12, 10)) if plot else None
        out.append(self.regress(X.replace('_AI', 'LH_PLUS_RH'), X, covariates, limits, plot=add_subplot(fh2, 2, 2, 1), cache_dir=cache_dir))
        out.append(self.regress(added_covariates_X[1], X, covariates, limits, plot=add_subplot(fh2, 2, 2, 2), cache_dir=cache_dir))

        out.append(self.regress(Y.replace('_AI', 'LH_PLUS_RH'), Y, covariates, limits, plot=add_subplot(fh2, 2, 2, 3), cache_dir=cache_dir))
        out.append(self.regress(added_covariates_Y[1], Y, covariates, limits, plot=add_subplot(fh2, 2, 2, 4), cache_dir=cache_dir))

        if np.all([o is None for o in out]):
            return None
        return out



def skip_key(key):
    return ('fuzzy' in key or
            key == 'MRI_cort_thick_ctx_mean_AI')


def skip_pairing(key1, key2):
    return (('AllFib' in key1 and 'AllFib' in key2) or
            ('AllFib' in key1 and 'DTI_' in key2) or
            ('AllFib' in key2 and 'DTI_' in key1) or
            np.any([t in key1 and t in key2 for t in ['SLF_', 'SCS_', '_Fx']]))


def find_one_relationship(all_data, key1, key2, covariates=[],
                          rsq_thresh=0., plot=False):
    print key1, key2, covariates

    # Limit to data without nan
    idx = np.ones(all_data[key1].shape, dtype=bool)
    for key in [key1, key2] + covariates:
        idx = np.logical_and(idx, np.logical_not(np.isnan(all_data[key])))
    if not np.any(idx):
        return None

    # Construct the data matrices
    data1 = all_data[key1][idx]
    data2 = all_data[key2][idx]
    X = [np.ones((idx.sum(),))]
    for key in [key2] + covariates:
        X.append(all_data[key][idx])

    # Do the regression, remove covariates.
    bestfit = np.linalg.lstsq(np.asarray(X).T, data1)[0]
    for ci, dat in enumerate(X[2:]):  # covariates
        data2 -= bestfit[ci + 2] * dat

    # Now, redo the regression on the residuals
    m, b, r, p, err = linregress(data1, data2)
    assert np.logical_not(np.isnan(r)), 'WTF, nan?'
    if r**2 < rsq_thresh:
        return None

    key = '%s vs. %s' % (key1, key2)

    # Small plot
    if plot:
        xlims = np.asarray([np.min(data1), np.max(data1)])
        ylims = np.asarray([np.min(data2), np.max(data2)])
        xlims += (np.diff(xlims) / 10. * np.asarray([-1, 1]))  # add space
        ylims += (np.diff(ylims) / 10. * np.asarray([-1, 1]))  # add space

        # Grab / create the axis handle.
        if isinstance(plot, bool):
            ax = plt.figure().gca()
        else:
            ax = plot

        # Plot axes
        ax.hold(True)
        if xlims[0] < 0 and xlims[1] > 0:
            ax.plot([0, 0], ylims, 'k--')
        if ylims[0] < 0 and ylims[1] > 0:
            ax.plot(xlims, [0, 0], 'k--')
        ax.scatter(data1, data2)
        ax.plot(xlims, b + m * xlims, 'r', linewidth=3.0)
        x_mean = data1.mean()
        ax.plot(x_mean, b + m * x_mean, 'g*', markersize=20.0)

        # Metadata
        ax.set_title("Significant at %.2e (r=%.3f)" % (p, r))
        ax.set_xlabel(key1)
        ax.set_ylabel(key2)
        ax.set_xlim(tuple(xlims))
        ax.set_ylim(tuple(ylims))

    return key, p, r


def search_one_asymmetry(**kwargs):
    sess = PINGDataSession()  # username/password stored in an env variable for me.
    sess.login()

    return sess.regress_multistep(X='MRI_cort_area_ctx_supramarginal_AI',
                                  Y='MRI_cort_area_ctx_rostralmiddlefrontal_AI',
                                  **kwargs)


def search_all_pairwise(plot=True, **kwargs):
    """For each pair of variables, look for a significant regression slope."""

    sess = PINGDataSession()  # username/password stored in an env variable for me.
    sess.login()

    all_data = compute_all_asymmetries(prefix=['MRI_cort_area', 'MRI_cort_thick',
                                               'MRI_subcort_vol', 'DTI_fiber_vol'])
    results = []

    keys = list(set(all_data.keys()) - set(('SubjID',)))
    for ki in range(len(keys)):
        key1 = keys[ki]
        if skip_key(key1):
            continue

        for ii in range(ki + 1, len(keys)):
            key2 = keys[ii]
            if skip_key(key2) or skip_pairing(key1, key2):
                continue

            try:
                result = sess.regress_multistep(key1, key2, plot=plot, **kwargs)
                results.append(result)
            except Exception as e:
                print "Exception: %s" % str(e)
            else:
                if plot:
                    if result is None:
                        plt.close()
                        plt.close()
                    else:
                        plt.show()


def search_all_vs_one(key, plot=False, **kwargs):
    """For each pair of variables, look for a significant regression slope."""

    sess = PINGDataSession()  # username/password stored in an env variable for me.
    sess.login()

    all_data = compute_all_asymmetries(prefix=['MRI_cort_area', 'MRI_cort_thick',
                                               'MRI_subcort_vol', 'DTI_fiber_vol'])
    results = []
    keys = list(set(all_data.keys()) - set(('SubjID',)))
    for loop_key in keys:
        results.append(sess.regress(key, loop_key, plot=plot, **kwargs))
        if plot:
            plt.show()

def search_all_vs_itself(covariates, plot=False, **kwargs):
    """For each pair of variables, look for a significant regression slope."""

    def add_subplot(fh, *args):
        return fh.add_subplot(*args) if plot else None

    sess = PINGDataSession()  # username/password stored in an env variable for me.
    sess.login()

    all_data = compute_all_asymmetries(prefix=['MRI_cort_area', 'MRI_cort_thick',
                                               'MRI_subcort_vol', 'DTI_fiber_vol'])
    results = []
    keys = list(set(all_data.keys()) - set(('SubjID',)))
    for X in keys:
        # Get relevant covariates
        try:
            added_covariates = sess.AI2flds(X)

            # Then regress on each side
            fh2 = plt.figure(figsize=(18, 6)) if plot else None
            sess.regress(X.replace('_AI', '_LH_PLUS_RH'),
                         X,
                         covariates=covariates, plot=add_subplot(fh2, 1, 3, 1),
                         **kwargs)
            sess.regress(added_covariates[1],
                         X,
                         covariates=covariates + ['MRI_cort_area_ctx_total_LH_PLUS_RH'],
                         plot=add_subplot(fh2, 1, 3, 2),
                         **kwargs)
            sess.regress('MRI_cort_area_ctx_total_LH_PLUS_RH',
                         X,
                         covariates=covariates,
                         plot=add_subplot(fh2, 1, 3, 3),
                         **kwargs)
        except Exception as e:
            print("Failed for %s (%s); moving on..." % (X, e))

        if plot:
            plt.show()