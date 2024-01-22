#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


process filter_testable_snps {
    conda params.conda
    tag "${chromosome}"
    label "med_mem"

    input:
        tuple val(chromosome), path(input_data)

    output:
        tuple val(chromosome), path(tested)

    script:
    tested = "${chromosome}.tested.bed"
    """
    
    python3 $moduleDir/bin/lrt.py \
        ${input_data} \
        ${params.samples_file} \
        ${chromosome} \
        --min_indivs_per_group ${params.min_indivs_per_group} \
        --min_groups ${params.min_groups} \
        --chrom ${chromosome} \
        --coverage_tr ${params.fdr_coverage_filter}
    """
}

process fit_random_effects_model {
    conda "/home/sabramov/miniconda3/envs/condR"
    tag "${chromosome}"
    label "med_mem"

    input:
        tuple val(chromosome), path(input_data)

    output:
        path name

    script:
    name = "${chromosome}.result.bed"
    """
    Rscript $moduleDir/bin/fitMixedEffectModel.R \
        ${input_data} \
        ${name}
    """
}

process differential_cavs {
    conda params.conda
    publishDir "${params.outdir}"
    label "high_mem"

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
            | filter_testable_snps
            | fit_random_effects_model

        pvals = out
            | collectFile(
                name: "${params.aggregation_key}.pvals.bed",
                sort: true,
                keepHeader: true,
                skip: 1
            )
        
        tested = filter_testable_snps.out 
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