#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"
params.pval_file = ""


process scan_with_moods {
    conda params.conda
    tag "${motif_id}"
    scratch true
    publishDir "${params.outdir}/moods_scans", pattern: "${name}"

    input:
        tuple val(motif_id), val(cluster_id), path(pwm_path)

    output:
        tuple val(motif_id), val(cluster_id), path(pwm_path), path(name)
    
    script:
    name = "${motif_id}.moods.log.bed.gz"
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
    | bgzip -c \
    > ${name}
    """
}


process motif_enrichment {
    publishDir "${params.outdir}/${params.aggregation_key}/${pval_file.simpleName}/counts", pattern: "${counts_file}"
    publishDir "${params.outdir}/${params.aggregation_key}/${pval_file.simpleName}/enrichment", pattern: "${enrichment_file}"
    scratch true
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), val(cluster_id), path(pwm_path), path(moods_file), path(pval_file)

    output:
        tuple val(motif_id), path(counts_file), path(pval_file), emit: counts
        tuple val(motif_id), path(enrichment_file), path(pval_file), emit: enrichment

    script:
    counts_file = "${motif_id}.counts.bed.gz"
    enrichment_file = "${motif_id}.enrichment.bed.gz"
    """
    zcat ${moods_file} | bedmap \
    --skip-unmapped \
    --sweep-all \
    --range 20 \
    --delim "|" \
    --multidelim ";" \
    --echo \
    --echo-map <(sort-bed ${pval_file}) \
     -    \
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
    tag "${motif_id}"
    conda params.conda
    scratch true

    input:
        tuple val(motif_id), path(counts_file), path(pval_file)

    output:
        tuple val(motif_id), path(name), path(pval_file)
    
    script:
    name = "${motif_id}.motif_stats.tsv"
    """
    python3 $moduleDir/bin/motif_stats.py ${motif_id} ${pval_file} ${counts_file} ${name}
    """


}

workflow calcEnrichment {
    take:
        args
    main:
        counts = motif_enrichment(args).counts //, motifs.map(it -> it[2]).collect())
        motif_ann = get_motif_stats(counts)
        .collectFile(
            storeDir: "${params.outdir}/${params.aggregation_key}/motif_stats",
            keepHeader: true, newLine: true, skip: 1) { it -> [[ "${item[2].simpleName}.bed", item[1].text]]}
    emit:
        counts
}

workflow motifEnrichment {
    take:
        pval_file
    main:
        motifs = Channel.fromPath(params.motifs_list)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.motif, row.cluster, file(row.motif_file)))
        moods_scans = scan_with_moods(motifs).combine(pval_file)
        enrichment = calcEnrichment(moods_scans)
    emit:
        enrichment
}

workflow {
    pvals = Channel.fromPath("${params.pval_file_dir}/*.bed")
        .map(it -> file(it))
    motifEnrichment(pvals)
}