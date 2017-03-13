#!/usr/bin/env Rscript

library(data.table)

unfiltered_results <- fread("output/unfiltered_results.tsv")

unfiltered_results[evalue < 1 & query_span > 500 & hit_span > 500, .(
    mito_hits = sum(grepl('mitochondria',
                          unique(hit_description),
                          ignore.case = TRUE))
), by = query_id][mito_hits > 0]


unfiltered_results[query_id == "scaffold333|size56947"]
