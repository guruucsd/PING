import datetime
import json
import requests

from ping import PINGSession


class GWASSession(PINGSession):

    def __init__(self, *args, **kwargs):
        super(GWASSession, self).__init__(*args, **kwargs)
        self.result_ids = None  # current dictionary of result IDs

    def get_results_ids(self, force=False):
        assert self.sess is not None, "Session must be started before calling get_results_ids."

        if self.result_ids is None or force:
            resp = self.sess.get('https://ping-dataportal.ucsd.edu/applications/GWAS/getListOfRuns.php?project_name=PING')
            self.result_ids = json.loads(resp.text)['runs']
            self.log('Fetched %d result ID(s)' % len(self.result_ids))
        return self.result_ids

    def get_results(self, **kwargs):
        assert self.sess is not None, "Session must be started before calling get_results."
        if 'measure' in kwargs:
            ids = [r['id'] for r in self.get_results_ids()
                           if r['yvalue'] == kwargs['measure']]
            return [self.get_results(id=id) for id in ids]
        url = 'https://ping-dataportal.ucsd.edu/applications/GWAS/downloadRunResult.php?project_name=PING&id=%s' % (kwargs['id'])
        resp = self.sess.get(url)
        if resp.text == u'':
            raise Exception('id not found: %s' % kwargs['id'])
        results = [lin.split('\t') for lin in resp.text.split('\n')]

        self.log('Fetched %d result(s) for id=%s.' % (len(results), kwargs['id']))
        return results

    def launch_run(self, measure, covariates=['Age_At_IMGExam']):
        assert self.sess is not None, "Session must be started before calling launch_run."
        try:
            self.log('Launching run for measure=%s ...' % measure)
            start_time = datetime.datetime.now()
            url = 'https://ping-dataportal.ucsd.edu/applications/GWAS/startRun.php?project_name=PING&com=%s&covariates=%s' % (
                measure, '+'.join(covariates))
            print url
            # resp = self.sess.get(url)
            time_diff = datetime.datetime.now() - start_time

            self.log('Completed run for measure=%s (time=%s)' % (measure, str(time_diff)))
            self.get_result_ids(force=True)
        except Exception as e:
            print e
        else:
            print resp.text

    def launch_and_retrieve_run(self, measure):
        self.launch_run(measure=measure)
        print self.get_results(measure=measure)

if __name__ == '__main__':
    sess = GWASSession()  # username/password stored in an env variable for me.
    sess.login()
    print sess.get_results_ids()
    print sess.get_results(measure=u'MRI_cort_area_ctx_frontalpole_AI')
    # print sess.launch_and_retrieve_run(measure='MRI_cort_area_ctx_rostralanteriorcingulate_AI')
