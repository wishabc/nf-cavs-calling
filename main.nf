#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom; addImputedCavs; aggregateAllPvalsNegbin; aggregateAllPvalsBinom } from "./pval_calc"

workflow {
    // Estimate BAD and call 1-st round CAVs
    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv()
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    // Exclude 1-st round CAVs 
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files)

    // Reestimate BAD, and add excluded SNVs
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_')
    new_badmap_join = filtered_vcfs.join(new_badmap)
    bad_intersections = intersectWithBadmap(new_badmap_join, 'nocavs_') 
    imputed_cavs = addImputedCavs(bad_intersections.join(intersect_files))

    // Collect statistics to fit negative binomial distriution
    merged_files = imputed_cavs
        .map(item -> item[1])
        .collectFile(name: 'badmaps.tsv',
            keepHeader: true,
            storeDir: "${params.outdir}/stats")   
    bad_merge_file = bads.combine(merged_files)
    weights_files = fitNegBinom(bad_merge_file)

    // Calculate binomial and negative binomial pvalues
    negbin_pvals = calcPvalNegbin(imputed_cavs, weights_files, 'nocavs_')[0]
    aggregateAllPvalsNegbin(negbin_pvals)
    binom_pvals = calcPvalBinom(imputed_cavs, 'nocavs_')[0]
    aggregateAllPvalsBinom(binom_pvals)
    
}