import datetime
import os
import simplejson
import time

from . import PINGSession


class GWASSession(PINGSession):

    def __init__(self, *args, **kwargs):
        super(GWASSession, self).__init__(*args, **kwargs)
        self.result_ids = None  # current dictionary of result IDs

    def get_results_ids(self, force=False):
        assert self.sess is not None, "Session must be started before calling get_results_ids."
        self.log("Retrieving all result IDs ...")

        if self.result_ids is None or force:
            resp_text = self.download_file('applications/GWAS/getListOfRuns.php?project_name={project_name}')
            self.result_ids = simplejson.loads(resp_text)['runs']
            self.log("Fetched %d result ID(s)" % len(self.result_ids))
        return self.result_ids

    def get_results(self, ids=None, measures=None, force=False, raw=False, output_dir='results/gwas'):
        """Can search by id or measure.
        Search by id returns a list of tuples,
            each tuple a triplet (SNP, effect size, pval)
        Search by measure returns a list of list of tuples,
            an outer list for each id corresponding to that measure.
        """
        assert int(ids is None) + int(measures is None) == 1, 'id or measure (and not both) must be passed.'
        assert self.sess is not None, "Session must be started before calling get_results."

        def fetch_and_cache_result(id):
            # Fetch remote file
            out_file = os.path.join(output_dir, '%s_GWAS.csv' % id)
            if os.path.exists(out_file) and not force:
                with open(out_file, 'r') as fp:
                    out_text = '\n'.join(fp.readlines())
            else:
                self.log("Retrieving results for id=%s ..." % id)
                out_text = self.download_file('applications/GWAS/downloadRunResult.php?project_name={project_name}&id=%s' % id,
                                              out_file=out_file)
                if out_text == '':
                    os.remove(out_file)
                    raise Exception('Results for GWAS ID=%s not found.' % id)
                else:
                    self.log("Wrote results to disk at %s." % out_file)

            # Parse the result
            results = [lin.split('\t') for lin in out_text.split('\n')]
            self.log("Fetched %d result(s) for id=%s." % (len(results), id))

            return out_text if raw else results


        # Map measure to (possibly multiple) ids, then recursively call.
        if measures is not None:
            ids = [r['id'] for r in self.get_results_ids()
                           if r['yvalue'] in measures]
        self.log("Retrieving results for measures %s (%d ids found) ..." % (
            measures, len(ids)))

        results = []
        for cur_id in ids:
            try:
                results.append(fetch_and_cache_result(id=cur_id))
            except Exception as e:
                print("Failed to get id=%s: %s" % (cur_id, str(e)))
                results.append(None)
        return results

    def launch_run(self, measure, covariates=['Age_At_IMGExam']):
        assert self.sess is not None, "Session must be started before calling launch_run."

        self.log("Launching run for measure=%s, covariates=%s ..." % (measure, str(covariates)))
        time.sleep(5.);
        try:
            # if 'y' != raw_input("Is this what you want to do? (y/N) => ").lower():
            #     return

            start_time = datetime.datetime.now()
            resp = self.make_request('applications/GWAS/startRun.php?project_name={project_name}&com=%s&covariates=%s' % (
                                         measure, '+'.join(covariates)))
            time_diff = datetime.datetime.now() - start_time

            self.log("Completed run for measure=%s (time=%s)" % (measure, str(time_diff)))
            self.log("Results info:", self.get_results_ids(force=True)[-1])
        except Exception as e:
            print(e)
        else:
            self.log(resp.text)

    def launch_and_retrieve_run(self, measure, covariates=['Age_At_IMGExam'], output_dir='results/gwas'):
        result_id = self.launch_run(measure=measure, covariates=covariates)
        results = self.get_results(measure=measure, force=True, raw=True, output_dir=output_dir)
        return results[-1]
