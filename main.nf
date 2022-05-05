#!/usr/bin/env nextflow
import { estimate_bad } from "./bad_estimation"
import { extract_and_filter } from "./extract_and_filter"
import { calc_pval } from "./calc_pval"

workflow {
    extract_and_filter() | estimate_bad | calc_pval
}