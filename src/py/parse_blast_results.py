#!/usr/bin/env python3

from Bio import SearchIO


###########################
# SEARCH/FILTER FUNCTIONS #
###########################

def mito_filter(Hit):
    return ('mitochondr' in Hit.description.lower() or
            'cytochrome oxidase' in Hit.description.lower())


def length_filter(hsp):
    return(hsp.query_span > 100)


def evalue_filter(hsp):
    return(hsp.evalue < 0.01)


def hsp_filters(hsp):
    return(length_filter(hsp) and evalue_filter(hsp))


def return_tsv_lines(hsp, scaffold_length):
    return('%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%f' %
           (hsp.query_id,         # 1
            hsp.hit_id,           # 2
            hsp.hit_description,  # 3
            scaffold_length,      # 4
            hsp.query_span,       # 5
            hsp.query_start,      # 6
            hsp.query_end,        # 7
            hsp.hit_span,         # 8
            hsp.hit_start,        # 9
            hsp.hit_end,          # 10
            hsp.evalue,           # 11
            hsp.bitscore))        # 12


# load results
blast_results = 'output/blast_contigs/blast_results.xml'
qresults = SearchIO.parse(blast_results, 'blast-xml')
qresults_list = [x for x in qresults]

# filter by e-value, query_span and text search
long_hsps = [x.hsp_filter(hsp_filters) for x in qresults_list]
mito_hits = [x.hit_filter(mito_filter) for x in long_hsps]
mito_candidates = [x for x in mito_hits if len(x) > 0]

# header for tsv output
tsv_header = '\t'.join(
    ['query_id',
     'hit_id',
     'hit_description',
     'scaffold_length',
     'query_span',
     'query_start',
     'query_end',
     'hit_span',
     'hit_start',
     'hit_end',
     'evalue',
     'bitscore'])

# output unfiltered results
tsv_lines = [tsv_header]
for qresult in qresults_list:
    scaffold_length = qresult.seq_len
    for hit in qresult:
        for hsp in hit:
            tsv_lines.append(return_tsv_lines(hsp, scaffold_length))

with open('output/blast_contigs/unfiltered_results.tsv', 'w') as f:
    f.write('\n'.join(tsv_lines))
    f.write('\n')


# output filtered results
tsv_lines = [tsv_header]
for qresult in mito_candidates:
    scaffold_length = qresult.seq_len
    for hit in qresult:
        for hsp in hit:
            tsv_lines.append(return_tsv_lines(hsp, scaffold_length))

with open('output/blast_contigs/filtered_results.tsv', 'w') as f:
    f.write('\n'.join(tsv_lines))
    f.write('\n')

