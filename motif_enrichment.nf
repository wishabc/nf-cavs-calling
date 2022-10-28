#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"
params.pval_file = ""


process scan_with_moods {
    conda params.conda
    tag "${motif_id}"
    scratch true

    input:
        tuple val(motif_id), val(cluster_id), val(pwm_path)

    output:
        tuple val(motif_id), val(cluster_id), path(name)
    
    script:
    name = "${motif_id}.moods.log"
    """
    { (cat ${params.bg_file} | head -n5 | tail -n +2 | cut -d" " -f2) || true; } > background_probs.py

    echo ${pwm_path}
    moods-dna.py --sep ";" -s ${params.alt_fasta_file} \
        --p-value ${params.pval_tr} --lo-bg `cat background_probs.py` \
        -m "${pwm_path}" -o moods.log
    
    cat moods.log | awk '{print \$1}' > chroms.txt

    cat moods.log \
    | cut -d";" -f2- \
    | sed 's/;\$//g' \
    | awk -v FS=";" -v OFS="\t" \
        '{ print \$2, \$2+length(\$5), \$1, \$4, \$3, \$5; }' \
    | sed 's/".pfm"/""/g' \
    | paste chroms.txt - \
    | sort-bed - \
    > ${name}
    """
}


process motif_enrichment {
    publishDir "${params.outdir}/motif_clean_figures", pattern: "${prefix}.pdf"
    publishDir "${params.outdir}/motif_enrichment", pattern: "${prefix}.txt"
    scratch true
    conda params.conda

    input:
        tuple val(motif_id), val(cluster_id), path(moods_file)
        val pval_file
        path all_pwms

    output:
        tuple val(motif_id), val(cluster_id), path(counts_file), path(enrichment_file)

    script:
    counts_file = "${motif_id}.counts.bed.gz"
    enrichment_file = "${motif_id}.enrichment.bed.gz"
    """
    bedmap \
    --skip-unmapped \
    --sweep-all \
    --range 20 \
    --delim "|" \
    --multidelim ";" \
    --echo \
    --echo-map <(zcat ${pval_file}) \
        ${moods_file} \
    | python $projectDir/bin/parse_variants_motifs.py \
        ${params.genome_fasta_file} \
        ./ \
    | sort-bed - \
    | bgzip -c \
    > ${counts_file}

    python3 ${projectDir}/bin/motif_enrichment.py  \
        ${pval_file} \
        ${counts_file} \
        ${motif_id} \
        ${cluster_id} | bgzip -c > ${enrichment_file}
    """
}

workflow motifEnrichment {
    take:
        pval_file
    main:
        motifs = Channel.fromPath(params.motifs_list)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.motif, row.cluster, row.motif_file))
        moods_scans = scan_with_moods(motifs)
        enrichment = motif_enrichment(moods_scans, pval_file, motifs.map(it -> it[2]).collect())
    emit:
        enrichment
}

workflow {
    motifEnrichment(params.pval_file)
}