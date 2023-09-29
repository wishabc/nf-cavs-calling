#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"

process calc_pval_binom {
    tag "${indiv_id}"
    conda params.conda

    input:
        tuple val(indiv_id), path(badmap_intersect_file)
        val prefix

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.pvalue.bed"
    """
    python3 $moduleDir/bin/calc_pval_binom.py \
        -I ${badmap_intersect_file} \
        -O ${name} \
        -a ${params.allele_tr} \
        --ct ${params.coverage_tr} \
        --recalc-w
    """
}


process aggregate_pvals {
    conda params.conda
    tag "${indiv_id}"

    input:
        tuple val(indiv_id), path(pval_file)
        val prefix

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.aggregation.bed"
    """
    python3 $moduleDir/bin/aggregation.py \
        -I ${pval_file} \
        -O ${name}
    """
}

process exclude_cavs {
    conda params.conda
    tag "${indiv_id}"
    
    input:
        tuple val(indiv_id), path(non_aggregated_snps)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.snps.bed"
    """
    cat ${non_aggregated_snps} 
        | awk -v OFS='\t' '((NR==1) || (\$NF <= 0.05)) {print}'
        | cut -f1-11 > ${name}
    """
}


process add_cavs {
    conda params.conda
    tag "${indiv_id}"
    scratch true

    input:
        tuple val(indiv_id), path(new_bad_annotated), path(old_badmap_annotated)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.added_cavs.intersect.bed"
    n_badmap = new_bad_annotated.name != 'empty' ? "-n ${new_badmap}" : ""
    """
    python3 $moduleDir/bin/add_cavs.py \
        -o ${old_badmap_annotated} \
        ${n_badmap} \
        --output not_sorted_cavs.bed

    head -1 not_sorted_cavs.bed > ${name}
    sort-bed not_sorted_cavs.bed >> ${name}
    """
}

workflow calcPvalBinom {
    take:
        data
        prefix
    main:
        pval_files = calc_pval_binom(data, prefix)
    emit:
        pval_files
}


workflow callCavsFirstRound {
    take:
        bad_annotations
        prefix
    main:
        no_cavs_snps = calcPvalBinom(bad_annotations, prefix)
            | exclude_cavs
    emit:
        no_cavs_snps
}

// DEFUNC, test workflow
workflow {
    extracted_vcfs = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> row.indiv_id)
        | unique()
        | map(indiv_id -> tuple(indiv_id, "${params.outdir}/snp_annotation/${indiv_id}*"))
        
    callCavsFirstRound(extracted_vcfs, "")
}
