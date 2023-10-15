#!/usr/bin/env nextflow

// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"


process random_sample {
    tag "seed:${step_start}-${step_start+params.samples_per_job}"
    conda params.conda
    label "sampling"

    input:
        tuple val(step_start), path(annotations_file), path(non_aggregated_file)

    output:
        path name

    script:
    name = "${step_start}.sampling.tsv"
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
        params.sampling_count = 100
        params.samples_per_job = 100
        total_count = params.sampling_count * params.samples_per_job 
        out = Channel.of(1..params.sampling_count) 
            | map(it -> it * params.samples_per_job)
            | combine(data)
            | random_sample
            | collectFile(
                storeDir: params.outdir,
                name: "subsampled.n${total_count}.tsv",
                keepHeader: true,
                skip: 1
            )
    emit:
        out
        
}


workflow {
    nonagr_files = Channel.fromPath(params.nonagr_pvals)
    params.cavs_annotation = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/snvs.annotations.bed.gz"
    Channel.fromPath(params.cavs_annotation)
        | combine(nonagr_files)
        | sampleVariants
}