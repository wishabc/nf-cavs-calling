#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
//include { extractAndFilter } from "./extract_and_filter"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom } from "./pval_calc"


workflow {
    filtered_vcfs_and_intersect = estimateBadByIndiv()
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_')
    new_badmap_join = filtered_vcfs.join(new_badmap)
    new_intersect_map = intersectWithBadmap(new_badmap_join, 'nocavs_')
    switch (params.strategy) {
        case 'binom':
            calcPvalBinom(new_intersect_map, 'nocavs_')
            break
        case 'negbin':
            weights_files = fitNegBinom(new_intersect_map)
            calcPvalNegbin(new_intersect_map, weights_files, 'nocavs_')
            break
        default:
            println("Provided strategy ${strategy} is not supported. Please use either 'binom' or 'negbin'")
    }
    
}