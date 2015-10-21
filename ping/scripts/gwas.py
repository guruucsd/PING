"""
Script for running GWAS app on PING data.
"""

import os
from collections import OrderedDict

import simplejson

from . import snps as snps_script
from ..ping.apps.gwas import GWASSession
from ..research.apps import ResearchArgParser


def get_chromosome_locations(snp_metadata):
    cur_chromosome = '0'
    chrom_locs = [0]
    for si, snp in enumerate(snp_metadata):
        if cur_chromosome == snp['chromosome']:
            continue
        cur_chromosome = snp['chromosome']
        chrom_locs.append(chrom_locs[-1] + int(snp_metadata[si - 2]['basepair']))
    chrom_locs.append(chrom_locs[-1] + int(snp_metadata[-1]['basepair']))

    return chrom_locs


def do_gwas(action, measures, covariates=None, output_format=None,
            data_dir='data', output_dir='data', force=False,
            username=None, passwd=None):
    if covariates is None:
        covariates = ['Age_At_IMGExam']

    sess = GWASSession(username=username, passwd=passwd, data_dir=data_dir)
    sess.login()

    # Get the data
    json_file = os.path.join(output_dir, 'GWAS_%s__%s.json' % (
        '_'.join(measures), '_'.join(covariates)))
    if action == 'launch':
        raw_data = [sess.launch_and_retrieve_run(measures=measures,
                                                 covariates=covariates)]
    elif force or not os.path.exists(json_file):
        raw_data = sess.get_results(measures=measures, force=True)

    else:
        raw_data = None

    # Get
    _, snp_metadata = snps_script.do_genes(action='view', gene='all', output_format='json', data_dir=data_dir, output_dir=output_dir)

    # Dump to json
    if raw_data is not None:
        data_dict = dict()
        for m, snps in zip(measures, raw_data):
            data_dict[m] = OrderedDict([(str(snp), dict(effect_size=float(str(es)),
                                                        pval=float(str(pval)),
                                                        metadata=snp_metadata.get(snp)))
                                        for snp, es, pval in snps])
        with open(json_file, 'w') as fp:
            simplejson.dump(data_dict, fp)

    # Display the data
    if output_format == 'json':
        pass  # already did the work

    elif output_format == 'flask':
        chrom_locations = get_chromosome_locations(snp_metadata.values())
        with open(os.path.join(output_dir, 'chrom_locs.json'), 'w') as fp:
            simplejson.dump(chrom_locations, fp)

        import flask
        app = flask.Flask(__name__)
        cur_dir = os.path.abspath(os.path.dirname(__file__))
        viz_dir = os.path.abspath(os.path.join(cur_dir, '..', 'viz', 'manhattan'))

        @app.route('/')
        def serve_default():
            return flask.send_from_directory(viz_dir, 'manhattan.html')

        @app.route('/data/<path:path>')
        def serve_data(path):
            return flask.send_from_directory(output_dir, path)

        @app.route('/<path:path>')
        def serve_everything(path):
            return flask.send_from_directory(viz_dir, path)

        app.run(port=5002)

    else:
        raise ValueError("Unknown output format: %s" % output_format)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Launch or view results of"
                               " a GWAS on the PING dataset.\n",
                               common_args=['username', 'passwd', 'force',
                                            'data-dir', 'output-dir'])
    parser.add_argument('action', choices=['display', 'launch'])
    parser.add_argument('measures', help="comma-separated list of "
                        "measures from the PING database,"
                        "including custom measures uploaded via upload.py")
    parser.add_argument('covariates', help="comma-separated list of "
                        "covariates from the PING database, "
                        "including custom measures uploaded via upload.py",
                        nargs='?', default='Age_At_IMGExam')
    parser.add_argument('--output-format', choices=['json', 'flask'],
                        nargs='?', default='json')
    args = parser.parse_args()
    args.measures = args.measures.split(',')
    args.covariates = args.covariates.split(',')

    do_gwas(**vars(args))
