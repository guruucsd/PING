"""
Access SNP information
"""
from ping.apps.snps import PINGSNPSession


def do_usage(args, error_msg=None):
    if error_msg is not None:
        print("*** ERROR *** : %s" % error_msg)
    print("\nUsage: %s {action} {SNP/gene}" % __file__)
    print("\tShow SNP=>gene or gene=>SNP mappings.")
    print("\n\taction: 'view' or 'download'")
    print("\t\tview: view mapping between SNP/gene, and any relevant GWAS results.")
    print("\t\tdownload: download all subject SNP data for the SNP/gene specified.")
    print("\t\t\tNOTE: Only 5000 total SNPs can be downloaded over the lifetime access to PING.")
    print("\tSNP/gene: case-sensitive text label; if it starts with 'rs', it is taken to be a SNP")


def do_snps(*args):

    if len(args) <= 1:
        do_usage(args, 'Not enough arguments.')

    elif args[0] not in ['view', 'download']:
        do_usage(args, 'Unknown command: %s' % args[0])

    elif args[1].startswith('rs'):
        # SNP => gene mapping
        snp = args[1]
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

        if args[0] == 'download':
            sess.login()
            sess.download_subject_snps([snp])

    else:
        # gene => SNP mapping
        gene_name = args[1]
        sess = PINGSNPSession()

        print("Loading genes...")
        gene_metadata = sess.get_gene_metadata(gene_name)
        print("Found %d genes for %s" % (len(gene_metadata), gene_name))
        print('\n'.join([str(gene) for gene in gene_metadata]))

        print("Loading snps...")
        all_snps = sess.get_snps_from_gene(gene_metadata)
        print("Found %d snps for %s" % (len(all_snps), gene_name))
        print('\n'.join([str(snp) for snp in all_snps]))

        if args[0] == 'download':
            sess.login()
            sess.download_subject_snps(all_snps)


if __name__ == '__main__':
    import sys
    do_snps(*sys.argv[1:])
