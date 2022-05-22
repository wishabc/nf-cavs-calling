#!/usr/bin/env nextflow
include { get_file_by_indiv_id } from "./helpers"
stats_dir = params.outdir + '/stats'

def get_snp_annotation_file_by_id(indiv_id) {
    return "${params.outdir}/snp_annotation/" + get_file_by_indiv_id(indiv_id, "intersect")
}

process collect_stats_for_negbin {

    publishDir stats_dir

    input:
        path bad_annotations
    output:
        path "${stats}"
    script:
    stats = './'
    """
    python3 $baseDir/bin/collect_nb_stats.py ${bad_annotations} ${stats}
    """
} 

process calculate_pvalue {

    tag "P-value calculation ${indiv_id}"
    publishDir params.outdir + "/pval_files"

    input:
        tuple val(indiv_id), path(badmap_intersect_file)
        path stats_file
        val strategy
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
    input:
        tuple val(indiv_id), path(pval_vcf)
        val strategy
        val output
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "aggregation-${strategy}")
    """
    python3 $baseDir/bin/aggregation.py -I ${pval_vcf} -O ${name}
    """
}

process exclude_cavs {
    input:
        tuple val(indiv_id), path(bad_annotations), path(agg_vcf)
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "filter")
    """
    python3 $baseDir/bin/filter_cavs.py -a ${agg_vcf} -b ${bad_annotations} -O ${name} --fdr ${param.excludeFdrTr}
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
    main:
        pval_files = calculate_pvalue(data, stats_file, 'negbin')
        agg_files = aggregate_pvals(pval_files, 'negbinom')
    emit:
        agg_files
}



workflow callCavsFromVcfs {
    take:
        bad_annotations
    main:
        // all_badmaps = bad_annotations
        //     .map{ it -> it[1] }
        //     .collectFile(name: 'bad_annotations_files.txt', newLine: true, storeDir: stats_dir)
        //stats_file = collect_stats_for_negbin(all_badmaps)
        agg_files = calcPvalBinom(bad_annotations, '')
        agg_file_cavs = bad_annotations.join(agg_files)
        no_cavs_snps = exclude_cavs(agg_file_cavs)
        //calcPvalNegbin(bad_annotations, stats_file)
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

workflow {
    callCavs()
}
