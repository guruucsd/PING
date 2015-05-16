"""
Allow access to generic SNP information, as well as
downloading subject SNP info from the SNP browser.
"""
import csv
import json
import os
import subprocess
import sys

import numpy as np
import six

from . import PINGSession

GENES = None
SNPS = None


class PINGSNPSession(PINGSession):
    """Allows a user to access generic PING data about gene/SNP relations,
    as well as download user-data for specific snps.
    """

    genes_metadata_file = 'csv/genes/PING_gene_annotate.json'
    SNP_metadata_file = 'csv/genes/PING_SNPs.txt'

    def get_genes_dict(self):
        global GENES
        if GENES is None:
            if not os.path.exists(self.genes_metadata_file):
                self.log("Downloading PING genes metadata...")
                self.download_file('data/PING/data_uncorrected/SNPs/PING_gene_annotate.json',
                                   out_file=self.genes_metadata_file)
            GENES = json.load(open(self.genes_metadata_file, 'r'))
            GENES = np.asarray(GENES['data'])
        return GENES

    def get_gene_metadata(self, gene_name):
        """Gene metadata contains gene names, chromosome location, and
        base pair range."""
        all_genes = self.get_genes_dict()

        idx = np.asarray([gene_name in g[8] for g in all_genes])
        return all_genes[idx]

    def get_snp_metadata(self, snp):
        if not os.path.exists(self.SNP_metadata_file):
            self.log("Downloading PING SNP metadata...")
            self.download_file('data/PING/data_uncorrected/SNPs/PING_SNPs.txt',
                               out_file=self.SNP_metadata_file)
        snp_reader = csv.reader(open(self.SNP_metadata_file, 'r'))
        next(snp_reader)  # skip header

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
        matched_genes = []
        for gene in all_genes:
            if (cur_chromosome == gene[1] and
                    int(gene[3]) - range/2 <= cur_basepair and
                    cur_basepair <= int(gene[4]) + range/2):
                matched_genes.append(gene)
        return matched_genes

    def get_snps_from_gene(self, gene, chromosome=None, range=100):
        """
        Given a gene name, or entry in the gene metadata dictionary,
        return all SNPs for the PING study.
        """

        # Prep the overall gene info
        if isinstance(gene, six.string_types):
            gene_metadata = self.get_gene_metadata(gene)
        else:
            gene_metadata = gene
        all_chromosomes = np.unique([g[1][3:] for g in gene_metadata])

        # Prep a stream of the SNP info
        if not os.path.exists(self.SNP_metadata_file):
            self.log("Downloading PING SNP metadata...")
            self.download_file('data/PING/data_uncorrected/SNPs/PING_SNPs.txt',
                               out_file=self.SNP_metadata_file)
        snp_reader = csv.reader(open(self.SNP_metadata_file, 'r'))
        next(snp_reader)  # skip header

        # Loop over each SNP in the file, saving matches
        matched_snps = []
        for row in snp_reader:
            # Format: ['SNP', 'Chromosome', 'Basepair', 'Allele1', 'Allele2']

            if row[1] not in all_chromosomes:
                continue  # Wrong chromosome, no need to process further!

            cur_chromosome = 'chr%s' % row[1]
            cur_basepair = int(row[2])

            # Search over the requested genes to see if we have a match.
            for gene in gene_metadata:
                if (cur_chromosome == (chromosome or gene[1]) and
                        int(gene[3]) - range/2 <= cur_basepair and
                        cur_basepair <= int(gene[4]) + range/2):
                    matched_snps.append(row)
                    break

        return matched_snps

    def download_subject_snps(self, all_snps):
        """Download actual data from subjects"""
        snp_ids = [s if isinstance(s, six.string_types) else s[0]
                   for s in all_snps]
        snp_txt = self.download_file('applications/SNPs/download.php?_v=&project_name={project_name}&snps=%s' % (
                                         '%0A'.join(snp_ids)),
                                     out_file='download/snps/%s.csv' % '_'.join(snp_ids))

        lines = snp_txt.split(' ')[4:]
        header = lines[0]

        return lines
