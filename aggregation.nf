#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
    | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
    | transpose()
}

process aggregate_pvals {
    conda params.conda
    tag "${indiv_id}"

    input:
        tuple val(indiv_id), path(pval_file)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.aggregation.bed"
    """
    python3 $moduleDir/bin/aggregation.py \
        -I ${pval_file} \
        -O ${name} \
        --max_coverage_tr ${params.fdr_coverage_filter} \
        --weights ${weights}
    """
}

process merge_files {
    conda params.conda
    scratch true
    tag "${group_key}"

    input:
        tuple val(group_key), path(files)

    output:
        tuple val(group_key), path(name)

    script:
    name = "${group_key}.sorted.bed"
    """
    echo "`head -n 1 ${files[0]}`\tgroup_id" > ${name}
    tail -n +2 -q ${files} | sed "s/\$/\t${group_key}/" | sort-bed - >> ${name}
    """
}

process pack_data {
    conda params.conda
    publishDir params.outdir
    scratch true

    input:
        path aggregated_variants
        path non_aggregated_variants
    
    output:
        path sorted_aggregated, emit: aggregated
        tuple path(sorted_non_aggregated), path("${sorted_non_aggregated}.tbi"), emit: non_aggregated
        path name, emit: stats
    
    script:
    sorted_aggregated = "aggregated.${params.aggregation_key}.bed"
    sorted_non_aggregated = "non_aggregated.${params.aggregation_key}.bed.gz"
    name = "cav_stats.${params.aggregation_key}.tsv"
    """
    head -1 ${aggregated_variants} > ${sorted_aggregated}
    sort-bed ${aggregated_variants} >> ${sorted_aggregated}

    head -1 ${non_aggregated_variants} > tmp.txt
    sort-bed ${non_aggregated_variants} >> tmp.txt
    bgzip -c tmp.txt > ${sorted_non_aggregated}
    tabix ${sorted_non_aggregated}

    python3 $moduleDir/bin/collect_stats.py \
        ${sorted_aggregated} \
        ${name} \
        --fdr ${params.fdr_tr}
    """

}


workflow aggregation {
    take:
        sample_split_pvals
    main:
        params.aggregation_key = params.aggregation_key ?: "all"
        if (params.aggregation_key != 'all') {
            sample_cl_correspondence = Channel.fromPath(params.samples_file)
                | splitCsv(header: true, sep: '\t')
                | map{ row ->
                    if (row.containsKey(params.aggregation_key)) {
                        return tuple(row.ag_id, row[params.aggregation_key])
                    } else {
                        throw new Exception("Column '${params.aggregation_key}' does not exist in the samples file '${params.samples_file}'")
                    }
                }
            pvals = sample_split_pvals
                | join(sample_cl_correspondence)
                | filter(it -> !it[2].isEmpty())
                | map(it -> tuple(it[2], it[1]))
                | set_key_for_group_tuple 
                | groupTuple() 
        } else {
            pvals = sample_split_pvals
                | map(it -> it[1])
                | collect(sort: true)
                | map(it -> tuple('all', it)) 
        }
        iter2_prefix = 'final'

        aggregated_merged = pvals
            | merge_files
            | aggregate_pvals
            | map(it -> it[1])
            | collectFile(
                skip: 1,
                keepHeader: true,
                name: "aggregated.bed",
            )
        non_aggregated_merged = merge_files.out
            | map(it -> it[1])
            | collectFile(
                skip: 1,
                keepHeader: true,
                name: "non_aggregated.bed",
            )
        packed = pack_data(aggregated_merged, non_aggregated_merged)
    emit:
        packed.aggregated
        packed.non_aggregated
}


workflow {
    by_sample_pvals = Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))

    aggregation(by_sample_pvals)
}
