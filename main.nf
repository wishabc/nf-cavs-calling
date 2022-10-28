#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; addImputedCavs } from "./pval_calc"
include { motifEnrichment } from "./motif_enrichment"


params.conda = "$moduleDir/environment.yml"


process sort_and_gzip {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path inp

    output:
        tuple path(name), path("${name}.tbi")

    script:
    name = "${inp.simpleName}.sorted.bed.gz"
    """
    sort-bed ${inp} | bgzip -c > ${name}
    tabix ${name}
    """
}


workflow test {
    binom_files = Channel.fromPath('/net/seq/data2/projects/sabramov/ENCODE4/cav-calling/babachi_1.5_common_final/output/final.pval_files_binom/*.pvalue.bed')
        .map(it -> file(it))
    binom_files.take(3).view()
    a = binom_files.collectFile(name: "all_variants.bed")
    a.view()
    sort_and_gzip(a) | motifEnrichment
}


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
    enrichment = sort_and_gzip(binom_pvals.collectFile(
        name: "all_variants.bed", skip:1
        )) | motifEnrichment
    
}