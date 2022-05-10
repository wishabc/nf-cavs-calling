#!/usr/bin/env nextflow
include { estimateBadAndIntersect } from "./bad_estimation"
include { extractAndFilter } from "./extract_and_filter"
//include { calcPval } from "./calc_pval"

workflow {
    extractAndFilter() | estimateBadAndIntersect
}