"""
Script for running GWAS app on PING data.
"""
import os

from ping.apps.gwas import GWASSession
from research.apps import ResearchArgParser


def do_gwas(action, measure, username=None, passwd=None):

    covariates = ['Age_At_IMGExam']  # ', 'DTI_fiber_vol_AllFibnoCC_AI', 'Gender', 'FDH_23_Handedness_Prtcpnt']  # MRI_cort_area_ctx_total_LH_PLUS_RH']
    print(action, measure, covariates)

    sess = GWASSession(username=username, passwd=passwd)
    sess.login()

    if action == 'display':
        print(sess.get_results(measure=measure, force=True))
    elif action == 'launch':
        print(sess.launch_and_retrieve_run(measure=measure, covariates=covariates))


if __name__ == '__main__':
    parser = ResearchArgParser(description="Launch or view results of"
                               " a GWAS on the PING dataset.\n",
                               common_args=['username', 'passwd'])
    parser.add_argument('action', choices=['display', 'launch'])
    parser.add_argument('measure', help="any measure from the PING database,"
                        "including custom measures uploaded via upload.py")
    args = parser.parse_args()
    do_gwas(**vars(args))
