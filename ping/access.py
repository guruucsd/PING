"""
Accessing PING data
"""
import md5
import os

import numpy as np
import pandas
import requests
import statsmodels.formula.api as smf

PING_DATA = None


class PINGSession(object):

    def __init__(self, username=None, passwd=None, verbosity=1):
        self.username = username or os.environ.get('PING_USERNAME')
        self.passwd = passwd or os.environ.get('PING_PASSWORD')

        if self.username is None:
            raise ValueError("username must be specified, or read from environment variables PING_USERNAME")
        elif '@' in self.username:
            raise ValueError('You must log in with your username, not email address, for these functions to work.')
        if self.passwd is None:
            raise ValueError("password must be specified, or read from environment variables PING_PASSWORD")

        self.sess = None  # http session with server
        self.result_ids = None  # current dictionary of result IDs
        self.verbosity = verbosity  # level of output

    def log(self, message, verbosity=1):
        if self.verbosity >= verbosity:
            print message

    def login(self):
        payload = {
            'username': self.username,
            'pw': md5.md5(self.passwd).hexdigest(),
            'ac': 'log',
            'url': ''}

        self.sess = requests.Session()
        resp = self.sess.post('https://ping-dataportal.ucsd.edu/applications/User/login.php',
                              data=payload)
        if 'Login to Data Portal' in resp.text or resp.url != 'https://ping-dataportal.ucsd.edu/index.php':
            self.sess = None
            raise Exception('Login failed.')
        else:
            self.log("Logged in as %s successfully." % self.username)

    def get_spreadsheet(self, out_file=None):
        url = 'https://ping-dataportal.ucsd.edu/applications/Documents/downloadDoc.php?project_name=PING&version=&file=../usercache_PING_%s.csv' % (
            self.username)
        resp = self.sess.get(url)
        out_text = str(resp.text)

        if out_file:
            with open(out_file, 'wb') as fp:
                fp.write(out_text)

        return out_text


def col2prop(col_name):
    return col_name.replace('-', '.')


def load_PING_data(scrub_fields=False, csv_path=None, username=None, passwd=None, force=False):
    global PING_DATA
    if PING_DATA is None or force:

        csv_path = csv_path or os.path.join('csv', 'PING_raw_data.csv')

        # Download data
        if not os.path.exists(csv_path):
            print("Downloading PING data...")
            sess = PINGSession(username=username, passwd=passwd)
            sess.login()
            sess.get_spreadsheet(out_file=csv_path)

        print("Loading PING data...")
        data = pandas.read_csv(csv_path)

        # Convert dots to underscores
        print("Converting PING data...")

        new_data = dict()
        for key in data.keys():
            if scrub_fields and '.' not in key:
                continue
            new_data[key.replace('.', '_')] = data[key].as_matrix()
        data = new_data

        # print("Regressing data on confounds...")
        # for key in data.keys():
        #     formula = ('%s ~ FDH_Highest_Education + FDH_3_Household_Income +'
        #                '     DeviceSerialNumber + GAF_africa + GAF_amerind +'
        #                '     GAF_eastAsia + GAF_oceania + GAF_centralAsia') % key.replace('.', '_')
        #     try:
        #         resid = smf.ols(formula, data=data).fit().resid
        #         data[key] = resid
        #     except Exception as e:
        #         print "Failed (%s): %s" % (key, e)
        PING_DATA = data
    return PING_DATA


def which_hemi(key):
    rh_markers = ['Right', '_rh_', '_R_']
    lh_markers = ['Left', '_lh_', '_L_']
    if np.any([m in key for m in rh_markers]):
        return 'rh'
    elif np.any([m in key for m in lh_markers]):
        return 'lh'
    else:
        return None


def get_twohemi_keys(prefix, all_keys):
    """Given a key prefix, get all keys that have left/right pairs."""

    # Massage inputs
    if not isinstance(prefix, list):
        prefix = [prefix]
    prefix = np.asarray(prefix)

    good_keys = []
    for prop_name in all_keys:
        if (which_hemi(prop_name) != 'rh' or
                not np.any([p in prop_name for p in prefix])):
            continue
        if 'vent' in prop_name.lower() and 'ventral' not in prop_name.lower():
            print("Skipping %s" % prop_name)
            continue
        good_keys.append(prop_name)
    return np.asarray(good_keys)



def get_bilateral_hemi_keys(data, prop_name):
    """ Given one hemisphere's property name,
    return both."""

    if 'Left' in prop_name or 'Right' in prop_name:
        left_prop_name = prop_name.replace('Right', 'Left')
        right_prop_name = prop_name.replace('Left', 'Right')
    elif '_lh_' in prop_name or '_rh_' in prop_name:
        left_prop_name = prop_name.replace('_rh_', '_lh_')
        right_prop_name = prop_name.replace('_lh_', '_rh_')
    elif '_L_' in prop_name or '_R_' in prop_name:
        left_prop_name = prop_name.replace('_R_', '_L_')
        right_prop_name = prop_name.replace('_L_', '_R_')
    else:
        raise ValueError("Unknown format for prop_name='%s'" % prop_name)

    return left_prop_name, right_prop_name
