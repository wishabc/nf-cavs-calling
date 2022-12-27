#!/usr/bin/env nextflow

stats_dir = "${params.outdir}/stats"
params.conda = "$moduleDir/environment.yml"

process calculate_pvalue {
    tag "${indiv_id}"
    publishDir "${params.outdir}/${output}pval_files_${strategy}"
    conda params.conda

    input:
        tuple val(indiv_id), path(badmap_intersect_file)
        path stats_file
        val strategy
        val output

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.pvalue.bed"
    """
    python3 $moduleDir/bin/calc_pval.py -I '${badmap_intersect_file}' \
         -O ${name} -s ${strategy} --stats-file ${stats_file}
    """
}


process aggregate_pvals {
    publishDir "${params.outdir}/${output}ag_files_${strategy}"
    conda params.conda
    tag "${indiv_id}"
    cpus 5

    input:
        tuple val(indiv_id), path(pval_vcf)
        val strategy
        val output

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.aggregation.bed"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    python3 $moduleDir/bin/aggregation.py -I '${pval_vcf}' \
        -O '${name}' --jobs ${task.cpus} \
        --ct ${params.fdr_cov_tr}
    """
}

process exclude_cavs {
    conda params.conda
    publishDir "${params.outdir}/excluded_cavs"
    tag "${indiv_id}"
    
    input:
        tuple val(indiv_id), path(bad_annotations), path(agg_vcf)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.snps.bed"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    python3 $moduleDir/bin/filter_cavs.py -a ${agg_vcf} \
     -b ${bad_annotations} -O ${name} \
     --fdr ${params.exclude_fdr_tr}
    """
}


process add_cavs {
    publishDir "${params.outdir}/added_cavs"
    conda params.conda
    tag "${indiv_id}"
    scratch true

    input:
        tuple val(indiv_id), path(new_badmap), path(old_badmap)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.added_cavs.intersect.bed"
    """
    python3 $moduleDir/bin/add_cavs.py -n ${new_badmap} \
     -o ${old_badmap} --output not_sorted_cavs.bed
    head -1 not_sorted_cavs.bed > ${name}
    sort-bed not_sorted_cavs.bed >> ${name}
    """
}

workflow addImputedCavs {
    take:
        data
    main:
        add_cavs(data)
    emit:
        add_cavs.out
}

workflow calcPvalBinom {
    take:
        data
        prefix
    main:
        pval_files = calculate_pvalue(data, "${projectDir}", 'binom', prefix)
        agg_files = aggregate_pvals(pval_files, 'binom', prefix)
    emit:
        pval_files
        agg_files
}


workflow callCavsFromVcfsBinom {
    take:
        bad_annotations
        prefix
    main:
        pval_agg_files = calcPvalBinom(bad_annotations, prefix)
        agg_files = pval_agg_files[1]
        agg_file_cavs = bad_annotations.join(agg_files)
        no_cavs_snps = exclude_cavs(agg_file_cavs)
    emit:
        no_cavs_snps
}

workflow callCavs {
    extracted_vcfs = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
        .map(row -> row.indiv_id)
        .unique()
        .map(indiv_id -> tuple(indiv_id, "${params.outdir}/snp_annotation/${indiv_id}*"))
        
    callCavsFromVcfsBinom(extracted_vcfs)
}


workflow {
    callCavs()
}
