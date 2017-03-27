#!/usr/bin/env python3

import tompytools
import tompltools
from Bio import SeqIO


# fasta length-finding function
def get_longest_scaffold(assembly_file):
    scaffold_length = sorted(
        len(rec) for rec in SeqIO.parse(assembly_file, 'fasta'))
    return(scaffold_length[-1])


def main():

    parsed_arguments = tompltools.parse_cli_arguments()

    assembly_files = tompytools.find_all(['noIUPAC.fasta'], 'output')

    # pick assembly with the longest scaffold
    longest_assembly = sorted(
        (get_longest_scaffold(x), x) for x in assembly_files)[-1]

    # get FASTA for the longest assembly
    fasta_file = longest_assembly[1]

    length_id = sorted(
        (len(rec), rec.id) for
        rec in SeqIO.parse(fasta_file, 'fasta'))

    longest_scaffold = length_id[-1][1]

    record_index = SeqIO.index(fasta_file, 'fasta')

    SeqIO.write(sequences=record_index[longest_scaffold],
                handle=parsed_arguments.output_fa,
                format='fasta')

if __name__ == '__main__':
    main()
