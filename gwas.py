"""
Script for running GWAS app on PING data.
"""
from ping.apps.gwas import GWASSession


def do_usage(args):
    print("\nUsage for %s:" % __file__)
    print("\t%s {action} {measure}" % __file__)
    print("\n\taction: 'display' or 'launch'")
    print("\t\tdisplay: show results from a previous run")
    print("\t\tlaunch: launch a new GWAS run")
    print("\tmeasure: any measure from the PING database, or measure uploaded via upload.py")
    print("\t\tThe measure will be regressed against variation in genes at each SNP.")
    print("\t\tAge_At_ImgExam will be used as the covariate.")
    print("\nExamples:")
    print("\t%s launch MRI_cort_area_ctx_total_LH_PLUS_RH" % __file__)
    print("\t\tThis launches a GWAS to search for genetic variation as a function of total cortical area.")
    print("\t%s display MRI_cort_area_ctx_total_LH_PLUS_RH" % __file__)
    print("\t\tThis downloads and displays the top 200 SNPs related to total cortical area.")


def do_gwas(*args):

    if len(args) != 2:
        do_usage(args)

    action = args[0]
    measure = args[1]
    covariates = ['Age_At_IMGExam']  # ', 'DTI_fiber_vol_AllFibnoCC_AI', 'Gender', 'FDH_23_Handedness_Prtcpnt']  # MRI_cort_area_ctx_total_LH_PLUS_RH']
    print(action, measure, covariates)

    sess = GWASSession()  # username/password stored in an env variable for me.
    sess.login()

    if action == 'display':
        print(sess.get_results(measure=measure, force=True))
    elif action == 'launch':
        print(sess.launch_and_retrieve_run(measure=measure, covariates=covariates))


if __name__ == '__main__':
    import sys
    do_gwas(*sys.argv[1:])
