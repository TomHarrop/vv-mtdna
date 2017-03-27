#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import tompltools
from Bio import Entrez
from Bio import SeqIO


def main():
    # require email address from environment
    if os.getenv('NCBI_EMAIL'):
        # identify user to Entrez
        Entrez.email = os.getenv('NCBI_EMAIL')
    else:
        raise EnvironmentError(
            'Email address not set.\n'
            'Specify email address with -e or --email.')

    # download COI gb
    output_fa = args.output_fa[0]
    vv_coi_gi = os.path.splitext(os.path.basename(output_fa))[0]
    handle = Entrez.efetch(
        db='nuccore',
        id=vv_coi_gi,
        rettype='gb',
        retmode='text')
    for rec in SeqIO.parse(handle, 'gb'):
        SeqIO.write(rec, 'data/' + rec.id + '.fasta', 'fasta')

if __name__ == '__main__':
    main()
