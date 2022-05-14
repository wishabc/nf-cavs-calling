#!/usr/bin/env nextflow
include { get_file_by_indiv_id } from "./helpers"
stats_dir = params.outdir + '/stats'

process collect_stats_for_negbin {
    publishDir stats_dir
    input:
        path bad_annotations
    output:
        path "stats"
    script:
    """
    echo collecting negbin stats
    """
}
// python3 collect_nb_stats.py ${bad_annotations} ${stats}

process calculate_pvalue {
    input:
        tuple val(indiv_id), path(badmap_intersect_file)
        path stats_file
        val strategy
    output:
        tuple val(indiv_id), path(name)

    script:
    name = get_file_by_indiv_id(indiv_id, "pvalue-${strategy}")
    """
    echo Calc bin pval ${indiv_id}
    """
}

process aggregate_pvals {
    input:
        tuple val(indiv_id), path(pval_vcf)
        val strategy
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "aggregation-${strategy}")
    """
    echo Aggregating pvals ${pval_vcf} ${strategy}
    """
}


workflow calcPvalBinom {
    take:
        data
    main:
        pval_files = calculate_pvalue(data, "", 'binom')
        agg_files = aggregate_pvals(pval_files, 'binom')
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
        all_badmaps = bad_annotations
            .map{ it -> it[1] }
            .collectFile(name: 'bad_annotations_files.txt', newLine: true, storeDir: stats_dir)
        stats_file = collect_stats_for_negbin(all_badmaps)
        calcPvalBinom(bad_annotations)
        calcPvalNegbin(bad_annotations, stats_file)
        
}


workflow callCavs {
    extracted_vcfs = Channel.fromPath(params.samplesFile)
        .splitCsv(header:true, sep:'\t')
        .map{ row -> tuple(row.indiv_id,
            "${params.outdir}/snp_annotation/${get_file_by_indiv_id(row.indiv_id, 'intersect')}") }
    
    callCavsFromVcfs(extracted_vcfs)
}

workflow {
    callCavs()
}
