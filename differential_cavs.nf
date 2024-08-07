#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


process filter_testable_snps {
    conda params.conda
    tag "${chromosome}"
    label "high_mem"

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

    // FIXME move to label
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
    scratch true

    input:
        path pvals
        path tested_snps
        val aggregation_key

    output:
        tuple path(pvals_new), path(tested_new), path(fit_fail)

    script:
    pvals_new = "differential_pvals.${aggregation_key}.bed"
    tested_new = "differential_tested.${aggregation_key}.bed"
    fit_fail = "differential_pvals.fit_fail.${aggregation_key}.bed"
    """
    python3 $moduleDir/bin/differential_cavs.py \
        ${tested_snps} \
        ${pvals} \
        tmp \
        --fdr ${params.diff_fdr_tr} \
        --coverage_tr ${params.fdr_coverage_filter}
    
    mv tmp.fit_fail.bed ${fit_fail}
    (head -1 tmp.pvals.bed && sort-bed tmp.pvals.bed | tail -n +2) > ${pvals_new}
    (head -1 tmp.tested.bed && sort-bed tmp.tested.bed | tail -n +2) > ${tested_new}
    """

}

workflow differentialCavs {
    take:
        data
    main:
        out = Channel.of(1..22)
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
        dif_cavs = differential_cavs(pvals, tested, params.aggregation_key)
    emit:
        dif_cavs
        tested
}

workflow {
    params.nonagr_pvals = "${params.outdir}/non_aggregated.${params.aggregation_key}.bed.gz"
    Channel.fromPath(params.nonagr_pvals)
        | differentialCavs
}