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
        -a ${params.allele_tr}
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
        tuple val(indiv_id), path(aggregated_snps), path(bad_annotations)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.snps.bed"
    """
    python3 $moduleDir/bin/filter_cavs.py \
        -a ${aggregated_snps} \
        -b ${bad_annotations} \
        -O ${name} \
        --fdr ${params.fdr_tr}
    """
}


process add_cavs {
    conda params.conda
    tag "${indiv_id}"
    scratch true

    input:
        tuple val(indiv_id), path(new_badmap), path(old_badmap)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.added_cavs.intersect.bed"
    n_badmap = new_badmap.name != 'empty' ? "-n ${new_badmap}" : ""
    """
    python3 $moduleDir/bin/add_cavs.py \
        -o ${old_badmap} \
        ${n_badmap} \
        --output not_sorted_cavs.bed

    head -1 not_sorted_cavs.bed > ${name}
    sort-bed not_sorted_cavs.bed >> ${name}
    """
}

workflow addExcludedCavs {
    take:
        data
    main:
        out = add_cavs(data)
    emit:
        out
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


workflow callCavsFromVcfsBinom {
    take:
        bad_annotations
        prefix
    main:
        pval_files = calcPvalBinom(bad_annotations, prefix)
        no_cavs_snps = aggregate_pvals(pval_files, prefix)
            | join(bad_annotations)
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
        
    callCavsFromVcfsBinom(extracted_vcfs, "")
}
