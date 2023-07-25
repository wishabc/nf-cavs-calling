#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


process differential_cavs {
    conda params.conda
    tag "${chromosome}"
    memory 100.GB

    input:
        tuple val(chromosome), path(input_data)

    output:
        tuple path(pvals), path(tested)

    script:
    pvals = "${chromosome}.pvals.bed"
    tested = "${chromosome}.tested.bed"
    """
    python3 $moduleDir/bin/anova.py \
        ${input_data} \
        ${chromosome} \
        --min_samples ${params.min_samples} \
        --min_groups ${params.min_groups} \
        --allele_tr ${params.allele_tr} \
        --chrom ${chromosome}
    """
}

workflow differentialCavs {
    take:
        data
    main:
        out = Channel.of(1..22, 'X', 'Y')
            | map(it -> "chr${it}")
            | combine(data)
            | differential_cavs

        pvals = out 
            | map(it -> it[0]) 
            | collectFile(
                name: "${params.aggregation_key}.pvals.bed",
                storeDir: params.outdir,
                sort: true,
                keepHeader: true,
                skip: 1
            )
        
        tested = out 
            | map(it -> it[1]) 
            | collectFile(
                name: "${params.aggregation_key}.tested.bed",
                storeDir: params.outdir,
                sort: true,
                keepHeader: true,
                skip: 1
            )
    emit:
        pvals
        tested
}

workflow {
    Channel.fromPath(params.nonagr_pvals)
        | differentialCavs
}