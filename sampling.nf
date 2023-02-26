#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"

process random_sample {
    tag "${step_start}"
    publishDir "${params.outdir}/sampling"
    conda params.conda

    input:
        val step_start

    output:
        path name

    script:
    name = "${step_start}.sampling.tsv"
    """
    python3 $moduleDir/bin/random_sample.py \
        -I ${params.sampling_file} \
        -O ${name} \
        -c ${params.context_file} \
        -m ${params.mutation_rates} \
        --start ${step_start} \
        --step ${params.samples_per_job}
    """
}

workflow {
    params.sampling_count = 1000
    params.samples_per_job = 10
    params.mutation_rates = "/home/sabramov/projects/mutationRates/output/mutation_rates/mut_rates.annotation.bed"
    params.context_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-annotations/context/output/variants_context.bed"
    params.sampling_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-cavs/output/pvals_nonaggregated.origin/Normal.sorted.bed"
    sampling_results = Channel.of(1..params.sampling_count - 1) 
        | map(it -> it * params.samples_per_job)
        | toInteger()
        | random_sample
}