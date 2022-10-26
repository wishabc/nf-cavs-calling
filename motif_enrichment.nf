#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"
params.pval_file = ""


process scan_with_sarus {
    conda params.conda
    tag "${motif_id}"

    input:
        tuple val(motif_id), path(pwm_path), path(precalc_thresholds)
        path seqs
    output:
        val(motif_id), path(sarus_log)
    
    script:
    sarus_log = "${motif_id}.sarus.log"
    """
    java -cp $moduleDir/sarus.jar ru.autosome.SARUS ${seqs} \
        ${pwm_path} -10000000 \
        --pvalues-file ${precalc_thresholds} \
        --threshold-mode score \
        --output-scoring-mode logpvalue > ${sarus_log}
    """
}


process motif_enrichment {
    publishDir "${params.outdir}/motif_clean_figures", pattern: "${prefix}.pdf"
    publishDir "${params.outdir}/motif_enrichment", pattern: "${prefix}.txt"

    conda params.conda

    input:
        tuple val(motif_id), val(motif)
        val pval_file

    output:
        tuple val(motif_id), path("${prefix}.pdf"), path("${prefix}.txt")

    script:
    prefix = "${motif_id}"
    """
    python3 ${projectDir}/bin/motif.py  \
	${pval_file}
	by_pfm/${prefix}.txt \
	${motif} \
	${prefix}
    """
}

workflow motifEnrichment {
    take:
        pval_file
    main:
        motifs = Channel.fromPath(params.motifs_list)
            .splitText()
            .map(it -> tuple(file(it).simpleName, it))
        enrichment = motif_enrichment(motifs, pval_file)
    emit:
        enrichment
}

workflow {
    motifEnrichment(params.pval_file)
}