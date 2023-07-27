#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"

process random_sample {
    tag "${step_start}"
    publishDir "${params.outdir}/sampling"
    conda params.conda

    input:
        val step_start
        tuple path(non_aggregated_file), path(annotations_file)

    output:
        path name

    script:
    name = "${step_start}.sampling.tsv"
    //cpg = params.noncpg ? '--noncpg' : ''
    """
    python3 $moduleDir/bin/random_sample.py \
        ${non_aggregated_file} \
        ${annotations_file} \
        ${name} \
        --start ${step_start} \
        --step ${params.samples_per_job}
    """
}

workflow sampleVariants {
    take:
        data
    main:
        params.sampling_count = 1000
        params.samples_per_job = 10
        steps = Channel.of(1..params.sampling_count) 
            | map(it -> it * params.samples_per_job)
        out = random_sample(steps, data)
    emit:
        out
        
}

workflow {
    params.noncpg = false
    
    params.sampling_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-cavs/output/pvals_nonaggregated.origin/Normal.sorted.bed"
    sampling_results = Channel.of(1..params.sampling_count) 
        | map(it -> it * params.samples_per_job)
        | toInteger()
        
}