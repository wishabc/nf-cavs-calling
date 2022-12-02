#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"
params.pval_file = ""


process scan_with_moods {
    conda params.conda
    tag "${motif_id}"

    input:
        tuple val(motif_id), val(cluster_id), path(pwm_path)

    output:
        tuple val(motif_id), val(cluster_id), path(pwm_path), path(name)
    
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
    publishDir "${params.outdir}/${pval_file.simpleName}/counts", pattern: "${counts_file}"
    publishDir "${params.outdir}/${pval_file.simpleName}/enrichment", pattern: "${enrichment_file}"
    scratch true
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), val(cluster_id), path(pwm_path), path(moods_file), path(pval_file)
        // path all_pwms

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
    --echo-map <(sort-bed ${pval_file}) \
        ${moods_file} \
    | python $projectDir/bin/parse_variants_motifs.py \
        ${params.genome_fasta_file} \
        ./ \
    | sort-bed - \
    | bgzip -c \
    > ${counts_file}
    if ! [ -f ${counts_file} ]; then
        exit 1
    fi

    python3 ${projectDir}/bin/motif_enrichment.py  \
        ${pval_file} \
        ${counts_file} \
        ${motif_id} \
        ${cluster_id} | bgzip -c > ${enrichment_file}
    """
}

process get_motif_stats {
    publishDir "${params.outdir}/${pval_file.simpleName}/motif_stats"
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), path(counts_file), path(pval_file)

    output:
        tuple val(motif_id), path(name)
    
    script:
    name = "${motif_id}.motif_stats.tsv"
    """
    python3 $moduleDir/bin/motif_stats.py ${motif_id} ${pval_file} ${counts_file} ${name}
    """


}

workflow motifEnrichment {
    take:
        pval_file
    main:
        motifs = Channel.fromPath(params.motifs_list)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.motif, row.cluster, file(row.motif_file)))
        moods_scans = scan_with_moods(motifs)
        args = moods_scans | combine(pval_file)
        enrichment = motif_enrichment(args) //, motifs.map(it -> it[2]).collect())
        arg = enrichment.map(it -> tuple(it[0], it[2])).combine(pval_file)
        motif_ann = get_motif_stats(arg)
            // .collectFile(name: "motif_stats.bed",
            //     storeDir: "${params.outdir}/${ag_key}",
            //     keepHeader: true, skip: 1)
    emit:
        enrichment
        motif_ann
}

workflow {
    motifEnrichment(Channel.of(params.pval_file))
}