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

from export_measures import compute_all_asymmetries
from ping import PINGSession
from search_for_relationships import find_one_relationship, print_legend, skip_key, skip_pairing


class PINGDataSession(PINGSession):

    def __init__(self, *args, **kwargs):
        super(PINGDataSession, self).__init__(*args, **kwargs)
        self.result_ids = None  # current dictionary of result IDs
        if '@' in self.username:
            raise Exception('You must log in with your username, not email address, for these functions to work.')

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
        out.append(self.regress(added_covariates_X[0], X, covariates, limits, plot=add_subplot(fh2, 2, 2, 1), cache_dir=cache_dir))
        out.append(self.regress(added_covariates_X[1], X, covariates, limits, plot=add_subplot(fh2, 2, 2, 2), cache_dir=cache_dir))

        out.append(self.regress(added_covariates_Y[0], Y, covariates, limits, plot=add_subplot(fh2, 2, 2, 3), cache_dir=cache_dir))
        out.append(self.regress(added_covariates_Y[1], Y, covariates, limits, plot=add_subplot(fh2, 2, 2, 4), cache_dir=cache_dir))

        if np.all([o is None for o in out]):
            return None
        return out


def search_one_asymmetry(**kwargs):
    sess = PINGDataSession()  # username/password stored in an env variable for me.
    sess.login()

    return sess.regress_multistep(X='MRI_cort_area_ctx_supramarginal_AI',
                                  Y='MRI_cort_area_ctx_rostralmiddlefrontal_AI',
                                  **kwargs)


def search_all_asymmetries(plot=True, **kwargs):
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

if __name__ == '__main__':
    try:
        plt.figure()
    except:
        plot=False
    else:
        plot=True
        plt.close()

    search_all_asymmetries(plot=plot)
    print_legend()

    if plot:
        plt.show()
