#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
//include { extractAndFilter } from "./extract_and_filter"
include { calcPvalForVcfs } from "./pval_calc"

workflow {
    intersect_map = estimateBadByIndiv()
    no_cavs_snps = calcPvalForVcfs(intersect_map) | excludeCavs
    new_badmap = estimateBad(no_cavs_snps)
    new_badmap_join = intersect_map.join(new_badmap)
    new_intersect_map = intersectWithBadmap(new_badmap_join)
    calcPvalForVcfs(new_intersect_map)
}