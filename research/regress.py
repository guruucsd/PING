"""
"""
from .computed_measures import compute_all_asymmetries
from ping.apps.regress import PINGDataSession


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
                print("Exception: %s" % str(e))  # force print
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