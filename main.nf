#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
//include { extractAndFilter } from "./extract_and_filter"
include { callCavsFromVcfs; calcPvalBinom } from "./pval_calc"

workflow {
    intersect_map = estimateBadByIndiv()
    no_cavs_snps = callCavsFromVcfs(intersect_map)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_badmap')
    new_badmap_join = intersect_map.join(new_badmap)
    new_intersect_map = intersectWithBadmap(new_badmap_join, 'nocavs_intersect')
    calcPvalBinom(new_intersect_map)
}