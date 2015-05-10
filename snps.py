"""
Access SNP information
"""
import sys

from ping.apps.snps import PINGSNPSession


# gwas
# similarity
# regress

"""
if __name__ == '__main__':
    if sys.argv[1] == 'parse':
        print(parse_PING_output('browser_result.txt' if len(sys.argv) < 3 else sys.argv[2]))
"""

sess = PINGSNPSession()

if sys.argv[1] in ['view', 'download']:
    if sys.argv[2].startswith('rs'):
        snp = sys.argv[2]

        print("Loading SNPS...")
        snp_metadata = sess.get_snp_metadata(snp)
        print("Found %d SNPs for %s" % (snp_metadata is not None, snp))
        print('\n'.join([str(snp) for snp in snp_metadata]))

        print("Loading genes...")
        all_genes = sess.get_genes_from_snp(snp)
        print("Found %d genes for %s" % (all_genes is not None, snp))
        print('\n'.join([str(gene) for gene in all_genes]))

        print("GWAS results:")
        print('\n'.join([str(gwas) for gwas in sess.snp_gwas_results(snp)]))

    else:
        gene_name = sys.argv[2]

        print("Loading genes...")
        gene_metadata = sess.get_gene_metadata(gene_name)
        print("Found %d genes for %s" % (len(gene_metadata), gene_name))
        print('\n'.join([str(gene) for gene in gene_metadata]))

        print("Loading snps...")
        all_snps = sess.get_snps_from_gene(gene_metadata)
        print("Found %d snps for %s" % (len(all_snps), gene_name))
        print('\n'.join([str(snp) for snp in all_snps]))

        if sys.argv[1] == 'download':
            sess.login()
            sess.download_snps(all_snps)

else:
    raise NotImplementedError('Unknown command: %s' % sys.argv[1])