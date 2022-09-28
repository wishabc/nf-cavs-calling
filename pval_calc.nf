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
    python3 $moduleDir/bin/calc_pval.py -I ${badmap_intersect_file} -O ${name} -s ${strategy} --stats-file ${stats_file} --es-method ${params.esMethod} ${params.recalcW ? "--recalc-w" : ""}
    """
}

process aggregate_pvals {
    publishDir "${params.outdir}/${output}ag_files_${strategy}"
    conda params.conda
    tag "${indiv_id}"
    cpus 3

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
    python3 $moduleDir/bin/aggregation.py -I ${pval_vcf} -O ${name} --jobs ${task.cpus} --mc ${params.fdrCovTr}
    """
}

process exclude_cavs {
    conda params.conda
    tag "${indiv_id}"
    
    input:
        tuple val(indiv_id), path(bad_annotations), path(agg_vcf)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.snps.bed"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    python3 $moduleDir/bin/filter_cavs.py -a ${agg_vcf} -b ${bad_annotations} -O ${name} --fdr ${params.excludeFdrTr}
    """
}

process fit_nb {
    publishDir "${params.outdir}/stats"
    tag "${bad}"
    cpus 2
    conda params.conda

    input:
        tuple val(bad), path(bad_annotations)

    output:
        path "BAD*/weights_*.tsv"

    script:
    out_path = './'
    """
    python3 $moduleDir/bin/collect_nb_stats.py -b ${bad_annotations} -O ${out_path} --bad ${bad}
    negbin_fit -O ${out_path} -m NB_AS -R 500 -r 500 --jobs ${task.cpus}
    """
}

process add_cavs {
    publishDir "${params.outdir}/added_cavs"
    conda params.conda
    tag "${indiv_id}"

    input:
        tuple val(indiv_id), path(new_badmap), path(old_badmap)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.added_cavs.intersect.bed"
    """
    python3 $moduleDir/bin/add_cavs.py -n ${new_badmap} -o ${old_badmap} --output ${name}
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
        pval_files = calculate_pvalue(data, "${projectDir}", 'binom', output)
        agg_files = aggregate_pvals(pval_files, 'binom', output)
    emit:
        pval_files
        agg_files
}

// workflow calcPvalNegbin {
//     take:
//         data
//         stats_file
//         output
//     main:
//         pval_files = calculate_pvalue(data, stats_file, 'negbin', output)
//         agg_files = aggregate_pvals(pval_files, 'negbin', output)
//     emit:
//         pval_files
// }



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

// workflow fitNegBinom {
//     take:
//         bad_merge_file
//     main:
//         fit_dir = fit_nb(bad_merge_file)
//         fit = fit_dir.collectFile(
//             name: 'negbin_fit_params.tsv',
//             keepHeader: true,
//             storeDir: stats_dir)
//     emit:   
//         fit.first()
// }

// workflow aggregateAllPvalsBeforeCavs {
//         take:
//         vcf_tuples
//     main:
//         all_pvals = Channel.of('ALL').combine(vcf_tuples.map(it -> it[1]).collectFile(
//             name: "ALL.pvals.binom.beforecavs.bed",
//             keepHeader: true,
//             storeDir: stats_dir).first())
//         aggregate_pvals(all_pvals, 'binom', 'beforecavs_')
// }
// workflow aggregateAllPvalsNegbin {
//     take:
//         vcf_tuples
//     main:
//         all_pvals = Channel.of('ALL').combine(vcf_tuples.map(it -> it[1]).collectFile(
//             name: "ALL.pvals.negbin.bed",
//             keepHeader: true,
//             storeDir: stats_dir).first())
//         aggregate_pvals(all_pvals, 'negbin', 'all_')
// }

// workflow aggregateAllPvalsBinom {
//     take:
//         vcf_tuples
//     main:
//         all_pvals = Channel.of('ALL').combine(vcf_tuples.map(it -> it[1]).collectFile(
//             name: "ALL.pvals.binom.bed",
//             keepHeader: true,
//             storeDir: stats_dir).first())
//         aggregate_pvals(all_pvals, 'binom', 'all_')
// }

workflow callCavs {
    extracted_vcfs = Channel.fromPath(params.samplesFile)
        .splitCsv(header:true, sep:'\t')
        .map(row -> row.indiv_id)
        .unique()
        .map(indiv_id -> tuple(indiv_id, "${params.outdir}/snp_annotation/${indiv_id}*"))
        
    callCavsFromVcfsBinom(extracted_vcfs)
}


workflow {
    callCavs()
}
