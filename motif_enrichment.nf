#!/usr/bin/env nextflow

include { filter_tested_variants } from "./main"
// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"


process calc_enrichment {
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), path(pval_file), path(counts_file)

    output:
        tuple val(motif_id), path(pval_file), path(name)
    
    script:
    name = "${motif_id}.stats.tsv"
    """
    python3 $moduleDir/bin/motif_stats.py  \
        ${pval_file} \
        ${counts_file} \
        ${name} \
        --flank_width ${params.flank_width} \
        --fdr ${params.fdr_tr}
    """
}


workflow cavsMotifEnrichment {
    take:
        data
    main:
        out = data 
            | motifCounts // motif_hits, motif_hits_index
            | combine(Channel.fromPath(params.result_file))
            | calc_enrichment
            | map(it -> it[2])
            | collectFile(
                storeDir: params.outdir,
                name: "motif_enrichmnent.${params.aggregation_key}.tsv",
                keepHeader: true,
                skip: 1
            )
    emit:
        out
}

workflow {
    params.motif_counts = "/net/seq/data2/projects/sabramov/ENCODE4/moods_scans.ref.0928/output/motif_counts"
    motif_counts = Channel.fromPath("${params.motif_counts}/*.counts.bed")
        | map(it -> tuple(it.name.replaceAll(".counts.bed", ""), it))

    data = Channel.fromPath(params.nonagr_pvals)
        | combine(motif_counts)
        | cavsMotifEnrichment

}
