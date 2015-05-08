"""
Script for running GWAS app on PING data.
"""
import sys

from ping.apps.gwas import GWASSession


action = 'display' if len(sys.argv) == 1 else sys.argv[1]
measure = 'MRI_cort_area_PC0' if len(sys.argv) <= 2 else sys.argv[2]
covariates = ['Age_At_IMGExam']  # ', 'DTI_fiber_vol_AllFibnoCC_AI', 'Gender', 'FDH_23_Handedness_Prtcpnt']  # MRI_cort_area_ctx_total_LH_PLUS_RH']
print action, measure, covariates

sess = GWASSession()  # username/password stored in an env variable for me.
sess.login()

if action == 'display':
    print sess.get_results(measure=measure, force=True)
elif action == 'launch':
    print sess.launch_and_retrieve_run(measure=measure, covariates=covariates)
