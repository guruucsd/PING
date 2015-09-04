"""
Access SNP information
"""
from argparse import ArgumentParser

from ping.apps.snps import PINGSNPSession


def do_snps(action, snp_gene):
    if snp_gene.startswith('rs'):
        # SNP => gene mapping
        snp = snp_gene
        sess = PINGSNPSession()

        print("Loading SNPS...")
        snp_metadata = sess.get_snp_metadata(snp)
        print("Found %d SNPs for %s" % (snp_metadata is not None, snp))
        print(snp_metadata)

        print("Loading genes...")
        all_genes = sess.get_genes_from_snp(snp)
        print("Found %d genes for %s" % (all_genes is not None, snp))
        print('\n'.join([str(gene) for gene in all_genes]))

        print("GWAS results:")
        print('\n'.join([str(gwas) for gwas in sess.snp_gwas_results(snp)]))

        if action == 'download':
            sess.login()
            sess.download_subject_snps([snp])

    else:
        # gene => SNP mapping
        gene_name = snp_gene
        sess = PINGSNPSession()

        print("Loading genes...")
        gene_metadata = sess.get_gene_metadata(gene_name)
        print("Found %d genes for %s" % (len(gene_metadata), gene_name))
        print('\n'.join([str(gene) for gene in gene_metadata]))

        print("Loading snps...")
        all_snps = sess.get_snps_from_gene(gene_metadata)
        print("Found %d snps for %s" % (len(all_snps), gene_name))
        print('\n'.join([str(snp) for snp in all_snps]))

        if action == 'download':
            sess.login()
            sess.download_subject_snps(all_snps)


if __name__ == '__main__':

    parser = ArgumentParser(description="Show SNP=>gene or gene=>SNP"
                                        " mappings.")
    parser.add_argument('action', choices=['view', 'download'])
    parser.add_argument('snp_gene', metavar="snp/gene",
                        help="case-sensitive text label; if it starts with"
                             " 'rs', it is taken to be a SNP")
    args = parser.parse_args()
    do_snps(**vars(args))
