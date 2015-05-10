"""
Accessing PING data
"""
import hashlib
import os

import numpy as np
import pandas
import requests
# import statsmodels.formula.api as smf

PING_DATA = None
raw_input = raw_input if 'raw_input' in dir() else input


class PINGSession(object):
    """
    """
    project_name = 'PING'
    base_url = 'https://ping-dataportal.ucsd.edu/'

    def __init__(self, username=None, passwd=None, verbosity=1):
        self.username = username or os.environ.get('PING_USERNAME')
        self.passwd = passwd or os.environ.get('PING_PASSWORD')

        if self.username is None:
            raise ValueError("username must be specified, or read from environment variables PING_USERNAME")
        elif '@' in self.username:
            raise ValueError('You must log in with your username, not email address, for these functions to work.')
        if self.passwd is None:
            raise ValueError("password must be specified, or read from environment variables PING_PASSWORD")

        self.sess = requests.Session()
        self.result_ids = None  # current dictionary of result IDs
        self.verbosity = verbosity  # level of output

    def log(self, message, verbosity=1):
        if self.verbosity >= verbosity:
            print(message)

    def make_url(self, rel_path):
        template_url = '{base_url}{rel_path}'.format(base_url=self.base_url,
                                                     rel_path=rel_path)
        return template_url.format(project_name=self.project_name)

    def make_request(self, rel_path, verb='get', **kwargs):
        url = self.make_url(rel_path)
        self.log("Downloading file from %s ..." % url)
        request_fn = getattr(self.sess, verb)
        resp = request_fn(url, **kwargs)
        self.log("Download completed.")

        return resp

    def download_file(self, rel_path, out_file=None, **kwargs):
        resp = self.make_request(rel_path=rel_path, **kwargs)
        out_text = resp.text

        if out_file:
            # Make the directory and dump the file.
            dir_path = os.path.dirname(out_file)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

            with open(out_file, 'w') as fp:
                fp.write(out_text)

        return out_text

    def login(self):
        payload = {
            'username': self.username,
            'pw': hashlib.md5(self.passwd.encode()).hexdigest(),
            'ac': 'log',
            'url': ''}

        self.log("Logging in as %s" % (self.username))
        resp = self.make_request(rel_path='applications/User/login.php',
                                 verb='post', data=payload)
        if 'Login to Data Portal' in resp.text or resp.url != self.make_url('index.php'):
            raise Exception('Login failed.')
        else:
            self.log("Logged in as %s successfully." % self.username)

    def expert_mode_script(self, X, Y, covariates=[]):
        """Produce the necessary 'expert' script"""
        return """dependent.measure     = "{Y}"
covariates.usr        = "{covariates}"
covariates.ses        = "FDH_Highest_Education + FDH_3_Household_Income"
covariates.dev        = "DeviceSerialNumber"
covariates.gaf        = "GAF_africa + GAF_amerind + GAF_eastAsia + GAF_oceania + GAF_centralAsia"
independent.variable  = "{X}"
smoothing.interaction = ""
""".format(X=X, Y=Y, covariates='+'.join(covariates))

    def regress(self, X, Y, covariates=[], **kwargs):
        """Do the regression remotely (via R), return the raw result.
        """
        cookie = np.random.randint(1000)
        payload = {
            '_v': '',
            'cookie': cookie,
            'user_name': self.username,
            'project_name': self.project_name,
            'command': '+'.join(covariates),
            'yvalue': Y,
            'functionOf': X,
            'interaction': '',
            'expert': self.expert_mode_script(X, Y, covariates)}
        payload.update(**kwargs)

        # Generate the regression and store the result.
        self.log("Computing regression for %s vs. %s..." % (X, Y))
        resp = self.make_request('applications/DataExploration/executeR.php',
                                 verb='post', data=payload)
        r_text = str(resp.text)
        return r_text

    def download_PING_spreadsheet(self, out_file=None):

        # Now access the PING data sheet
        out_text = self.download_file(
            rel_path='applications/Documents/downloadDoc.php?project_name={project_name}&version=&file=../usercache_PING_%s.csv' % (
                self.username),
            out_file=out_file)

        return out_text


def get_lh_prop_name(prop_name):
    return prop_name.replace('_rh_', '_lh_').replace('_Right_', '_Left_').replace('_R_', '_L_')


def get_nonhemi_prop_name(prop_name):
    return prop_name.replace('_rh_', '_').replace('_Right_', '_').replace('_R_', '_')


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
            sess.download_PING_spreadsheet(out_file=csv_path)

        print("Loading PING data...")
        try:
            data = pandas.read_csv(csv_path, low_memory=False)
        except ValueError as ve:
            # Corrupt spreadsheet. Re-GET
            print("Error loading the PING data: %s" % ve)
            yn = raw_input("The PING spreadsheet is corrupt. Delete and download? (y/N) > ")
            if yn.lower() == 'y':
                os.remove(csv_path)
                return load_PING_data(scrub_fields=scrub_fields, csv_path=csv_path,
                                      username=username, passwd=passwd, force=True)
            raise ve

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
        #         print("Failed (%s): %s" % (key, e))
        PING_DATA = data

    return PING_DATA


def get_tbx_data(data=None):
    if data is None:
        data = load_PING_data()
    tbx_data = dict()
    for key in data.keys():
        if not key.startswith('TBX_'):
            continue
        if data[key].dtype.name in ['string', 'object']:
            continue
        tbx_data[key] = data[key][good_subj_idx]
    return tbx_data


def get_fdh_data(data=None):
    if data is None:
        data = load_PING_data()
    fdh_data = dict()
    for key in data.keys():
        if not key.startswith('FDH_'):
            continue
        if data[key].dtype.name in ['string', 'object']:
            continue
        fdh_data[key] = data[key][good_subj_idx]
    return fdh_data

def which_hemi(key):
    rh_markers = ['Right', '_rh_', '_R_']
    lh_markers = ['Left', '_lh_', '_L_']
    if np.any([m in key for m in rh_markers]):
        return 'rh'
    elif np.any([m in key for m in lh_markers]):
        return 'lh'
    else:
        return None


def get_twohemi_keys(all_keys, prefix=None):
    """Given a key prefix, get all keys that have left/right pairs."""

    # Massage inputs
    if prefix is not None:
        if not isinstance(prefix, list):
            prefix = [prefix]
        prefix = np.asarray(prefix)

    good_keys = []
    for prop_name in all_keys:
        if which_hemi(prop_name) != 'rh':
            continue
        elif prefix is not None and not np.any([p in prop_name for p in prefix]):
            continue
        if 'vent' in prop_name.lower() and 'ventral' not in prop_name.lower():
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
