#!/usr/bin/env nextflow
include { get_file_by_indiv_id; get_stats_dir } from "./helpers"


stats_dir = get_stats_dir()

def get_snp_annotation_file_by_id(indiv_id) {
    return "${params.outdir}/snp_annotation/" + get_file_by_indiv_id(indiv_id, "intersect")
}

process calculate_pvalue {

    tag "P-value calculation ${indiv_id}"
    publishDir "${params.outdir}/${output}pval_files"

    input:
        tuple val(indiv_id), path(badmap_intersect_file)
        path stats_file
        val strategy
        val output
    output:
        tuple val(indiv_id), path(name)

    script:
    name = get_file_by_indiv_id(indiv_id, "pvalue-${strategy}")
    """
    python3 $baseDir/bin/calc_pval.py -I ${badmap_intersect_file} -O ${name} -s ${strategy} --stats-file ${stats_file}
    """
}

process aggregate_pvals {
    publishDir "${params.outdir}/${output}ag_files"

    cpus 2
    input:
        tuple val(indiv_id), path(pval_vcf)
        val strategy
        val output
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "aggregation-${strategy}")
    """
    python3 $baseDir/bin/aggregation.py -I ${pval_vcf} -O ${name} --jobs ${task.cpus} --mc ${params.fdrCovTr}
    """
}

process exclude_cavs {
    publishDir "${params.outdir}/excluded_cavs"
    input:
        tuple val(indiv_id), path(bad_annotations), path(agg_vcf)
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "filter")
    """
    python3 $baseDir/bin/filter_cavs.py -a ${agg_vcf} -b ${bad_annotations} -O ${name} --fdr ${params.excludeFdrTr}
    """
}

process fit_nb {
    publishDir stats_dir
    conda "/home/sabramov/miniconda3/envs/negbinfit"
    input:
        tuple val(bad), path(bad_annotations)
    output:
        path "BAD*/weights_*.tsv"
    script:
    out_path = './'
    """
    python3 $baseDir/bin/collect_nb_stats.py -b ${bad_annotations} -O ${out_path} --bad ${bad}
    negbin_fit -O ${out_path} -m NB_AS
    """
}


workflow calcPvalBinom {
    take:
        data
        output
    main:
        pval_files = calculate_pvalue(data, params.outdir, 'binom', output)
        agg_files = aggregate_pvals(pval_files, 'binom', output)
    emit:
        agg_files
}

workflow calcPvalNegbin {
    take:
        data
        stats_file
        output
    main:
        stats_list = data.combine(stats_file).map(it -> it[1])
        pval_files = calculate_pvalue(data, stats_list, 'negbin', output)
        agg_files = aggregate_pvals(pval_files, 'negbinom', output)
    emit:
        agg_files
}



workflow callCavsFromVcfsBinom {
    take:
        bad_annotations
    main:
        agg_files = calcPvalBinom(bad_annotations, '')
        agg_file_cavs = bad_annotations.join(agg_files)
        no_cavs_snps = exclude_cavs(agg_file_cavs)
    emit:
        no_cavs_snps
}

workflow callCavs {
    extracted_vcfs = Channel.fromPath(params.samplesFile)
        .splitCsv(header:true, sep:'\t')
        .map(row -> row.indiv_id)
        .distinct()
        .map( indiv_id -> tuple(indiv_id, get_snp_annotation_file_by_id(indiv_id)))
        
    callCavsFromVcfs(extracted_vcfs)
}

workflow fitNegBinom {
    take:
        bad_merge_file
    main:
        fit_dir = fit_nb(bad_merge_file)
        fit = fit_dir.collectFile(
            name: 'negbin_fit_params.tsv',
            keepHeader: true,
            storeDir: stats_dir)
    emit:   
        fit
}



workflow {
    callCavs()
}
