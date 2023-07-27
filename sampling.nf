#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"

workflow {
    params.noncpg = false
    
    params.sampling_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-cavs/output/pvals_nonaggregated.origin/Normal.sorted.bed"
    sampling_results = Channel.of(1..params.sampling_count) 
        | map(it -> it * params.samples_per_job)
        | toInteger()
        
}