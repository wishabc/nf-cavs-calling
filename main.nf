#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom } from "./pval_calc"


workflow {
    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv()
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_')
    new_badmap_join = filtered_vcfs.join(new_badmap)
    bad_intersections = intersectWithBadmap(new_badmap_join, 'nocavs_') 

    merged_files = bad_intersections
        .map(item -> item[1])
        .collectFile(name: 'badmaps.tsv',
            keepHeader: true,
            storeDir: stats_dir)   
    bad_merge_file = bads.combine(merged_files).first()
    bad_merge_file.view()

    weights_files = fitNegBinom(bad_merge_file)
    //calcPvalNegbin(new_intersect_map, weights_files, 'nocavs_')
    calcPvalBinom(new_intersect_map, 'nocavs_')
    
}