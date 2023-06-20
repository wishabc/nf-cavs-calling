#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process differential_cavs {
    conda params.conda
    publishDir "${params.outdir}/differential"

    input:
        path input_data

    output:
        tuple path(anova), path(melt)

    script:
    anova = "${params.aggregation_key}.differential_pvals.bed"
    melt = "${params.aggregation_key}.differential_tested.bed"
    """
    python3 $moduleDir/bin/anova.py \
        ${input_data} \
        ${params.aggregation_key} \
        --ct ${params.coverage_tr} \
        --min_samples ${params.min_samples} \
        --min_groups ${params.min_groups} \
        
    """
}

workflow {
    differential_cavs(file(params.nonagr_pvals))
}