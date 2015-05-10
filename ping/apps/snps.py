"""
Allow access to generic SNP information, as well as
downloading subject SNP info from the SNP browser.
"""
import csv
import json
import subprocess
import sys

import numpy as np
import six

from ..access import PINGSession

GENES = None
SNPS = None


class PINGSNPSession(PINGSession):
    """Allows a user to access generic PING data about gene/SNP relations,
    as well as download user-data for specific snps.
    """

    def get_genes_dict(self):
        global GENES
        if GENES is None:
            GENES = json.load(open('csv/genes/PING_gene_annotate.json'))
            GENES = np.asarray(GENES['data'])
        return GENES

    def get_gene_metadata(self, gene_name):
        """Gene metadata contains gene names, chromosome location, and
        base pair range."""
        all_genes = self.get_genes_dict()

        idx = np.asarray([gene_name in g[8] for g in all_genes])
        return all_genes[idx]

    def get_snp_metadata(self, snp):
        snp_reader = csv.reader(open('csv/genes/PING_SNPs.txt'))
        snp_reader.next()  # skip header

        for row in snp_reader:
            # Format: ['SNP', 'Chromosome', 'Basepair', 'Allele1', 'Allele2']
            if row[0] == snp:
                return row
        return None

    def snp_gwas_results(self, snp):
        import glob
        import os
        import re

        matches = []
        for csv_file in glob.glob('download/gwas/*.csv'):
            with open(csv_file, 'r') as fp:
                for line in fp:
                    if re.search(snp, line):
                        matches.append('%s: %s' % (csv_file, line[:-1]))
        return matches

    def get_genes_from_snp(self, snp, range=100):
        snp_metadata = self.get_snp_metadata(snp)
        if snp_metadata is None:
            raise Exception("Could not find snp %s" % snp)

        cur_basepair = int(snp_metadata[2])
        cur_chromosome = 'chr%s' % snp_metadata[1]
        all_genes = self.get_genes_dict()
        for gene in all_genes:
            if (cur_chromosome == gene[1] and
                    int(gene[3]) - range/2 <= cur_basepair and
                    cur_basepair <= int(gene[4]) + range/2):
                return gene
        return None

    def get_snps_from_gene(self, gene, chromosome=None, range=100):
        """
        Given a gene name, or entry in the gene metadata dictionary,
        return all SNPs for the PING study.
        """
        if isinstance(gene, six.string_types):
            gene_metadata = self.get_gene_metadata(gene)
        else:
            gene_metadata = gene

        matched_snps = []
        snp_reader = csv.reader(open('csv/genes/PING_SNPs.txt'))
        snp_reader.next()  # skip header

        all_chromosomes = [int(g[1][3:]) for g in gene_metadata]
        min_chromosome = np.min(all_chromosomes)
        max_chromosome = np.max(all_chromosomes)

        for row in snp_reader:
            # Format: ['SNP', 'Chromosome', 'Basepair', 'Allele1', 'Allele2']

            if int(row[1]) < min_chromosome or max_chromosome < int(row[1]):
                continue

            cur_chromosome = 'chr%s' % row[1]
            cur_basepair = int(row[2])

            for gene in gene_metadata:
                if (cur_chromosome == (chromosome or gene[1]) and
                        int(gene[3]) - range/2 <= cur_basepair and
                        cur_basepair <= int(gene[4]) + range/2):
                    matched_snps.append(row)
                    print("Hit #%4d (chromosome %s)" % (
                        len(matched_snps), cur_chromosome))
                    break
        print(matched_snps)
        return matched_snps

    def download_snps(self, all_snps):
        """Download actual data from subjects"""
        template_url = self.make_url('applications/SNPs/download.php?_v=&project_name={project_name}&snps=%s')
        fetch_url = template_url % ('%0A'.join([s[0] for s in all_snps]))
        print(fetch_url)
        import pdb; pdb.set_trace()

    def parse_PING_output(self, output_file):
        """DEPRECATED: Use download_snps to download and parse.
        Parse data from the """

        with open(output_file) as fp:
            alltext = fp.readlines()[0]
        lines = alltext.split(' ')[4:]
        header = lines[0]
        import pdb; pdb.set_trace()
        return lines

