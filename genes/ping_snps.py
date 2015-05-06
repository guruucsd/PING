import csv
import json
import numpy as np
import sys


GENES = None

def get_gene_metadata(gene_name):
    global GENES
    if GENES is None:
        GENES = json.load(open('PING_gene_annotate.json'))
        # ["#hg19.knownGene.name", "hg19.knownGene.chrom", "hg19.knownGene.strand",
        #  "hg19.knownGene.txStart", "hg19.knownGene.txEnd", "hg19.knownGene.cdsStart",
        #  "hg19.knownGene.cdsEnd", "hg19.kgXref.geneSymbol", "hg19.kgXref.description"]
        GENES = np.asarray(GENES['data'])

    idx = np.asarray([gene_name in g[8] for g in GENES])
    return GENES[idx]


def get_snps(gene_metadata, chromosome=None, range=100):
    matched_snps = []
    snp_reader = csv.reader(open('PING_SNPs.txt'))
    snp_reader.next()  # skip header

    all_chromosomes = [int(g[1][3:]) for g in gene_metadata]
    min_chromosome = np.min(all_chromosomes)
    max_chromosome = np.max(all_chromosomes)

    for row in snp_reader:
        #['SNP', 'Chromosome', 'Basepair', 'Allele1', 'Allele2']

        if int(row[1]) < min_chromosome or max_chromosome < int(row[1]):
            continue

        cur_chromosome = 'chr%s' % row[1]
        cur_basepair = int(row[2])

        for gene in gene_metadata:
            #print cur_chromosome, gene[1], gene[3], cur_basepair, gene[4]
            if (cur_chromosome == (chromosome or gene[1]) and
                 int(gene[3]) - range/2 <= cur_basepair and
                 cur_basepair <= int(gene[4]) + range/2):
                matched_snps.append(row)
                print "hit #%4d (chromosome %s)" % (len(matched_snps), cur_chromosome)
                break
    print matched_snps
    return matched_snps


def download_snps(all_snps):
    template_url = 'https://ping-dataportal.ucsd.edu/applications/SNPs/download.php?_v=&project_name=PING&snps=%s'
    fetch_url = template_url % ('%0A'.join([s[0] for s in all_snps]))
    print fetch_url


def parse_PING_output(output_file):
    with open(output_file) as fp:
        alltext = fp.readlines()[0]
    lines = alltext.split(' ')[4:]
    header = lines[0]
    import pdb; pdb.set_trace()
    return lines


if __name__ == '__main__':
    if sys.argv[1] == 'parse':
        print parse_PING_output('browser_result.txt' if len(sys.argv) < 3 else sys.argv[2])

    elif sys.argv[1] in ['view', 'download']:
        gene_name = sys.argv[2]

        print("Loading genes...")
        gene_metadata = get_gene_metadata(gene_name)
        print("Found %d genes for %s" % (len(gene_metadata), gene_name))
        print gene_metadata

        print("Loading snps...")
        all_snps = get_snps(gene_metadata)
        print("Found %d snps for %s" % (len(all_snps), gene_name))

        if sys.argv[1] == 'download':
            download_snps(all_snps)
