"""
Access SNP information
"""

import copy
import os.path

import simplejson

from ..ping.apps.snps import PINGSNPSession
from ..research.apps import ResearchArgParser


def do_snps(action, snp, output_format='print',
            output_dir='.', username=None, passwd=None):
    sess = PINGSNPSession(username=username, passwd=passwd)

    print("Loading SNPS...")
    snp_metadata = sess.get_snp_metadata(snp)
    print("Found %d SNPs for %s" % (snp_metadata is not None, snp))
    if output_format == 'print':
        print(snp_metadata)
    else:
        raise NotImplementedError()

    print("Loading genes...")
    all_genes = sess.get_genes_from_snp(snp)
    print("Found %d genes for %s" % (all_genes is not None, snp))
    if output_format == 'print':
        print('\n'.join([str(gene) for gene in all_genes]))
    else:
        raise NotImplementedError()

    print("GWAS results:")
    if output_format == 'print':
        print('\n'.join([str(gwas) for gwas in sess.snp_gwas_results(snp)]))
    else:
        raise NotImplementedError()

    if action == 'download':
        sess.login()
        sess.download_subject_snps([snp])

    return snp_metadata, all_genes


def do_genes(action, gene, output_format='print',
             output_dir='.', username=None, passwd=None):

    # gene => SNP mapping
    sess = PINGSNPSession(username=username, passwd=passwd)

    # Find matching genes
    print("Loading genes...")
    gene_metadata = sess.get_gene_metadata(gene if gene != 'all' else None)
    print("Found %d genes for %s" % (len(gene_metadata), gene))

    # Output matching genes
    gene_json = os.path.join(output_dir, 'genes_%s.json' % gene)
    if output_format == 'print':
        print('\n'.join([str(gene) for gene in gene_metadata]))
    elif output_format == 'json':
        with open(gene_json, 'w') as fp:
            simplejson.dump(gene_metadata, fp)

    # Find matching SNPs
    print("Loading snps...")
    if gene == 'all':
        snp_metadata = sess.get_snps()
    else:
        snp_metadata = sess.get_snps_from_gene(gene_metadata)

    # Output matching SNPs
    snps_json = os.path.join(output_dir, 'SNPs_%s.json' % gene)
    print("Found %d snps for %s" % (len(snp_metadata), gene))
    if output_format == 'print':
        print('\n'.join([str(snp) for snp in snp_metadata]))
    elif output_format == 'json':
        with open(snps_json, 'w') as fp:
            simplejson.dump(snp_metadata, fp)

    # Download subject data
    if action == 'download':
        sess.login()
        sess.download_subject_snps(snp_metadata)

    return gene_metadata, snp_metadata


def do_snps_datadump(action, snp_gene, output_format='print',
                     output_dir='.', username=None, passwd=None):

    if snp_gene.startswith('rs'):
        do_snps(action=action, snp=snp_gene, output_format=output_format,
                output_dir=output_dir, username=username, passwd=passwd)

    else:
        do_genes(action=action, gene=snp_gene, output_format=output_format,
                 output_dir=output_dir, username=username, passwd=passwd)


if __name__ == '__main__':

    parser = ResearchArgParser(description="Show SNP=>gene or gene=>SNP mappings.",
                               common_args=['username', 'passwd', 'output-dir'])
    parser.add_argument('action', choices=['view', 'download'])
    parser.add_argument('snp_gene', metavar="snp/gene",
                        help="case-sensitive text label; if it starts with"
                             " 'rs', it is taken to be a SNP")
    parser.add_argument('--output-format', choices=['print', 'json'],
                        nargs='?', default='print')
    args = parser.parse_args()
    do_snps_datadump(**vars(args))
