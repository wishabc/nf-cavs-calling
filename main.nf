#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom; addImputedCavs; aggregateAllPvalsNegbin; aggregateAllPvalsBinom } from "./pval_calc"

workflow {
    // Estimate BAD and call 1-st round CAVs
    iter1_prefix = 'iter1.'

    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv(iter1_prefix)
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files, iter1_prefix)


    iter2_prefix = 'final.'
    // Reestimate BAD, and add excluded SNVs
    iter2_intersections = estimateBad(no_cavs_snps, iter2_prefix)
    imputed_cavs = addImputedCavs(iter2_intersections.join(intersect_files))
    binom_pvals = calcPvalBinom(imputed_cavs, iter2_prefix)
    // Collect statistics to fit negative binomial distriution
    // merged_files = imputed_cavs
    //     .map(item -> item[1])
    //     .collectFile(name: 'badmaps.tsv',
    //         keepHeader: true,
    //         storeDir: "${params.outdir}/stats")   
    // bad_merge_file = bads.combine(merged_files)
    // weights_files = fitNegBinom(bad_merge_file)

    // Calculate binomial and negative binomial pvalues
    // negbin_pvals = calcPvalNegbin(imputed_cavs, weights_files, 'nocavs_')[0]
    // aggregateAllPvalsNegbin(negbin_pvals)
    
    //aggregateAllPvalsBinom(binom_pvals)
    
}