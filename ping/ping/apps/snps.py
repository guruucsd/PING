"""
Allow access to generic SNP information, as well as
downloading subject SNP info from the SNP browser.
"""
import csv
import os
import simplejson
import subprocess
import sys
from collections import OrderedDict

import numpy as np
import six

from . import PINGSession

GENES = None
SNPS = None


class PINGSNPSession(PINGSession):
    """Allows a user to access generic PING data about gene/SNP relations,
    as well as download user-data for specific snps.
    """

    genes_metadata_file = 'genes/PING_gene_annotate.json'
    SNP_metadata_file = 'genes/PING_SNPs.txt'
    SNP_HEADER = ['name', 'chromosome', 'basepair', 'allele1', 'allele2']
    GENE_HEADER = ['name', 'chromosome', 'strand', 'txStart', 'txEnd',
                   'cdsStart', 'cdsEnd', 'geneSymbol', 'description']

    def __init__(self, *args, **kwargs):
        super(PINGSNPSession, self).__init__(*args, **kwargs)
        self.genes_metadata_file = os.path.join(self.data_dir, PINGSNPSession.genes_metadata_file)
        self.SNP_metadata_file = os.path.join(self.data_dir, PINGSNPSession.SNP_metadata_file)

    def get_genes_dict(self):
        global GENES
        if GENES is None:
            if not os.path.exists(self.genes_metadata_file):
                self.log("Downloading PING genes metadata to %s..." % self.data_dir)
                self.download_file(os.path.join(self.data_dir, 'PING/data_uncorrected/SNPs/PING_gene_annotate.json'),
                                   out_file=self.genes_metadata_file)
            GENES = simplejson.load(open(self.genes_metadata_file, 'r'))
            GENES = np.asarray(GENES['data'])
            GENES = [dict(zip(self.GENE_HEADER, G)) for G in GENES]
        return GENES

    def get_gene_metadata(self, gene_name=None):
        """Gene metadata contains gene names, chromosome location, and
        base pair range."""
        all_genes = self.get_genes_dict()

        if gene_name is not None:
            return [g for g in all_genes if gene_name == g['geneSymbol']]
        else:
            return all_genes

    def _get_snp_stream(self):
        # Prep a stream of the SNP info
        if not os.path.exists(self.SNP_metadata_file):
            self.log("Downloading PING SNP metadata...")
            self.download_file(os.path.join(self.data_dir, 'PING/data_uncorrected/SNPs/PING_SNPs.txt'),
                               out_file=self.SNP_metadata_file)
        snp_reader = csv.reader(open(self.SNP_metadata_file, 'r'))
        next(snp_reader)  # skip header

        return snp_reader

    def get_snp_metadata(self, snp=None):
        for row in self._get_snp_stream():
            # Format:
            if snp is row[0] == snp:
                return dict(zip(self.SNP_HEADER, row))
        return None

    def snp_gwas_results(self, snp):
        import glob
        import os
        import re

        matches = []
        for csv_file in glob.glob(os.path.join(self.data_dir, 'gwas/*.csv')):
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
            if (cur_chromosome == gene['chromosome'] and
                    int(gene['txStart']) - range/2 <= cur_basepair and
                    cur_basepair <= int(gene['txEnd']) + range/2):
                matched_genes.append(gene)
        return matched_genes

    def get_snps(self):
        # Loop over each SNP in the file, saving matches
        matched_snps = OrderedDict()   # name, chromosome, basepair, allele1, allele2
        for row in self._get_snp_stream():
            snp_vals = dict(zip(self.SNP_HEADER, row))
            matched_snps[snp_vals['name']] = snp_vals
        return matched_snps

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
        all_chromosomes = np.unique([g['chromosome'][3:] for g in gene_metadata])

        # Loop over each SNP in the file, saving matches
        matched_snps = OrderedDict()
        for row in self._get_snp_stream():
            # Format: ['name', 'Chromosome', 'Basepair', 'Allele1', 'Allele2']
            if row[1] not in all_chromosomes:
                continue  # Wrong chromosome, no need to process further!
            snp_vals = dict(zip(self.SNP_HEADER, row))
            cur_chromosome = 'chr%s' % snp_vals['chromosome']
            cur_basepair = int(snp_vals['basepair'])

            # Search over the requested genes to see if we have a match.
            for gene in gene_metadata:
                if (cur_chromosome == (chromosome or gene['chromosome']) and
                        int(gene['txStart']) - range/2 <= cur_basepair and
                        cur_basepair <= int(gene['txEnd']) + range/2):
                    matched_snps[snp_vals['name']] = snp_vals
                    break

        return matched_snps

    def download_subject_snps(self, all_snps):
        """Download actual data from subjects"""
        snp_ids = [s if isinstance(s, six.string_types) else s[0]
                   for s in all_snps]
        snp_txt = self.download_file('applications/SNPs/download.php?_v=&project_name={project_name}&snps=%s' % (
                                         '%0A'.join(snp_ids)),
                                     out_file=os.path.join(self.data_dir, 'snps', '%s.csv' % '_'.join(snp_ids)))

        lines = snp_txt.split(' ')[4:]
        header = lines[0]

        return lines
