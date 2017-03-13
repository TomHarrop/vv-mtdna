#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import tompltools
from Bio import Entrez
from Bio import SeqIO


def main():
    # require email address from cli
    args = tompltools.parse_cli_arguments()
    if not args.email:
        raise ValueError('Email address is required')

    # identify user to Entrez
    Entrez.email = args.email

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
