#!/usr/bin/env nextflow
// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"


process calc_enrichment {
    tag "${motif_id}"
    conda params.conda
    label "med_mem"

    input:
        tuple val(motif_id), path(counts_file), path(pval_file)

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
            | calc_enrichment
            | map(it -> it[2])
            | collectFile(
                storeDir: params.outdir,
                name: "motif_enrichment.${params.aggregation_key}.tsv",
                keepHeader: true,
                skip: 1
            )
    emit:
        out
}

workflow {
    params.agr_pvals = "${params.outdir}/aggregated.${params.aggregation_key}.bed"
    motif_counts = Channel.fromPath("${params.motif_counts_ref}/*.counts.bed")
        | map(it -> tuple(it.name.replaceAll(".counts.bed", ""), it, file(params.agr_pvals)))
        | cavsMotifEnrichment

}
