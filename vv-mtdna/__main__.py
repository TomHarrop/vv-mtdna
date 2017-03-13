#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import ruffus
import tompltools
import os


#############
# FUNCTIONS #
#############
def find_all(names, path):
    path_contents = [
        (dirpath, filenames) for dirpath, dirnames, filenames
        in os.walk('data', followlinks=True)]
    all_path_files = []
    for dirpath, filenames in path_contents:
        for filename in filenames:
            all_path_files.append(os.path.join(dirpath, filename))
    names_path_matches = []
    for file in all_path_files:
        if any([os.path.basename(file) in x for x in names]):
            names_path_matches.append(file)
    return names_path_matches


############
# PIPELINE #
############
def main():

    # parse CLI
    parser = ruffus.cmdline.get_argparse(
        description='Vv mtDNA assembly pipeline')
    parser.add_argument('--email', '-e',
                        help='Email address, reported to NCBI',
                        type=str,
                        dest='email')

    options = parser.parse_args()

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines['main']

    # TEST FUNCTION
    test_job_function = tompltools.generate_job_function(
        job_script='src/sh/io_parser',
        job_name='test',
        verbose=True)

    # download COI seed file
    download_coi_fasta = main_pipeline.originate(
        name='download_coi_fasta.py',
        task_func=tompltools.generate_job_function(
            job_type='originate',
            job_script='src/py/download_coi_fasta.py',
            job_name='download_coi_fasta.py',
            extras=True),
        output='data/GU207861.1.fasta',
        extras=['-e ' + options.email])

    # define files
    sample_list = 'data/samples.txt'

    with open(sample_list, 'r') as f:
        csvreader = csv.reader(f)
        next(csvreader)
        file_list = {x[0]: [x[1], x[2]] for x in csvreader}

    pe_filenames = file_list['pe']
    mp_filenames = file_list['mp']

    pe_files = find_all(pe_filenames, 'data')

    # filter out weird hidden directories, what are these anyway?
    pe_files_filtered = [x for x in pe_files if '/.' not in x]

    # load files into ruffus
    raw_fq_files = main_pipeline.originate(
        name='raw_fq_files',
        task_func=os.path.isfile,
        output=pe_files_filtered)

    # trim adaptors
    trim_bbduk = main_pipeline.merge(
        name='trim_bbduk',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/trim_bbduk',
            job_name='trim_bbduk',
            cpus_per_task=8),
        input=raw_fq_files,
        output='output/trim_bbduk/pe_trimmed.fastq.gz')

    # run mitobim
    main_pipeline.transform(
        name='run_mitobim',
        # task_func=tompltools.generate_job_function(
        #     job_script='src/sh/run_mitobim',
        #     job_name='run_mitobim'),
        task_func=test_job_function,
        input=trim_bbduk,
        add_inputs=ruffus.add_inputs(download_coi_fasta),
        filter=ruffus.formatter(),
        output='output/mitobim/mitobim.log.txt')

    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph(
        "ruffus/flowchart.pdf", "pdf",
        pipeline_name="Vv mtDNA assembly pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=32)

if __name__ == '__main__':
    main()
