#!/usr/bin/env nextflow
include { estimate_bad } from "./bad_estimation"
include { extract_and_filter } from "./extract_and_filter"
include { calc_pval } from "./calc_pval"

workflow {
    extract_and_filter() | estimate_bad | calc_pval
}