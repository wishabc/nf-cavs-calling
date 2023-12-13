#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


process anova {
    conda params.conda
    tag "${chromosome}"
    label "high_mem"

    input:
        tuple val(chromosome), path(input_data)

    output:
        tuple path(pvals), path(tested)

    script:
    pvals = "${chromosome}.pvals.bed"
    tested = "${chromosome}.tested.bed"
    """
    
    python3 $moduleDir/bin/lrt.py \
        ${input_data} \
        ${chromosome} \
        --min_samples ${params.min_samples} \
        --min_groups ${params.min_groups} \
        --chrom ${chromosome} \
        --coverage_tr ${params.fdr_coverage_filter}
    """
}

process differential_cavs {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path pvals
        path tested_snps

    output:
        tuple path(pvals_new), path(tested_new)

    script:
    pvals_new = "differential_pvals.${params.aggregation_key}.bed"
    tested_new = "differential_tested.${params.aggregation_key}.bed"
    """
    python3 $moduleDir/bin/differential_cavs.py \
        ${tested_snps} \
        ${pvals} \
        tmp.bed \
        --fdr ${params.diff_fdr_tr} \
        --es ${params.es_tr}

    head -1 tmp.bed > ${pvals_new}
    sort-bed tmp.bed >> ${pvals_new}

    head -1 ${tested_snps} > ${tested_new}
    sort-bed ${tested_snps} >> ${tested_new}
    """

}

workflow differentialCavs {
    take:
        data
    main:
        out = Channel.of(1..22, 'X', 'Y')
            | map(it -> "chr${it}")
            | combine(data)
            | anova

        pvals = out 
            | map(it -> it[0]) 
            | collectFile(
                name: "${params.aggregation_key}.pvals.bed",
                sort: true,
                keepHeader: true,
                skip: 1
            )
        
        tested = out 
            | map(it -> it[1]) 
            | collectFile(
                name: "${params.aggregation_key}.tested.bed",
                sort: true,
                keepHeader: true,
                skip: 1
            )
        dif_cavs = differential_cavs(pvals, tested)
    emit:
        dif_cavs
        tested
}

workflow {
    Channel.fromPath(params.nonagr_pvals)
        | differentialCavs
}