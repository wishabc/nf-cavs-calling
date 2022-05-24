#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
//include { extractAndFilter } from "./extract_and_filter"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom } from "./pval_calc"


workflow {
    intersect_map = estimateBadByIndiv()
    no_cavs_snps = callCavsFromVcfsBinom(intersect_map)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_')
    new_badmap_join = intersect_map.join(new_badmap)
    new_intersect_map = intersectWithBadmap(new_badmap_join, 'nocavs_')
    switch (params.strategy) {
        case 'binom':
            calcPvalBinom(new_intersect_map, 'nocavs_')
            break
        case 'negbin':
            merged_files = new_intersect_map
            .map(item -> item[1])
            .collectFile(name: 'badmaps.tsv', keepHeader: true, cache: false)
            weights_files = fitNegBinom(merged_files)
            calcPvalNegbin(new_intersect_map, weights_files, 'nocavs_')
            break
        default:
            println("Strategy ${strategy} provided is not supported. Please use either 'binom' or 'negbin'")
    }
    
}