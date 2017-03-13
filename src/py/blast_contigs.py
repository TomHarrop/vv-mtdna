#!/usr/bin/env python3

import os
from Bio.Blast.Applications import NcbiblastnCommandline

scaffolds_file = 'data/scaffolds_sorted.fasta'
blast_results = 'output/blast_contigs/blast_results.xml'


def main():

    # make sure BLASTDB env is set
    if os.getenv('BLASTDB'):
        print('Using BLAST database folder: %s' % os.getenv('BLASTDB'))
    else:
        raise EnvironmentError("BLASTDB environment variable not set")

    # run blast
    blastn_cline = NcbiblastnCommandline(
        query=scaffolds_file,
        out=blast_results,
        outfmt=5,
        db='nt',
        evalue=1,
        max_target_seqs=10,
        num_threads=8,
        task='megablast')
    print(blastn_cline)
    blastn_cline()

if __name__ == '__main__':
    main()
