#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom; addImputedCavs; aggregateAllPvalsNegbin; aggregateAllPvalsBinom } from "./pval_calc"
include { get_stats_dir } from "./helpers"

stats_dir = get_stats_dir()

workflow {
    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv()
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_')
    new_badmap_join = filtered_vcfs.join(new_badmap)
    bad_intersections = intersectWithBadmap(new_badmap_join, 'nocavs_') 
    imputed_cavs = addImputedCavs(bad_intersections.join(intersect_files))

    merged_files = imputed_cavs
        .map(item -> item[1])
        .collectFile(name: 'badmaps.tsv',
            keepHeader: true,
            storeDir: stats_dir)   
    bad_merge_file = bads.combine(merged_files)
    weights_files = fitNegBinom(bad_merge_file)

    negbin_pvals = calcPvalNegbin(imputed_cavs, weights_files, 'nocavs_')
    //aggregateAllPvalsNegbin(negbin_pvals.map(it -> it[0]))
    binom_pvals = calcPvalBinom(imputed_cavs, 'nocavs_')
    //aggregateAllPvalsBinom(binom_pvals.map(it -> it[0]))
    
}