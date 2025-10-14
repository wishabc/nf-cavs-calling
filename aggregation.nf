#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
    | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
    | transpose()
}

process aggregate_pvals {
    conda params.conda
    tag "${group_key}"

    input:
        tuple val(group_key), path(pval_file)

    output:
        tuple val(group_key), path(name)

    script:
    name = "aggregation.bed"
    """
    python3 $moduleDir/bin/aggregation.py \
        -I ${pval_file} \
        -O ${name} \
        --max_coverage_tr ${params.fdr_coverage_filter} \
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
    name = "merged.sorted.bed"
    """
    echo "`head -n 1 ${files[0]}`\tgroup_id" > ${name}
    tail -n +2 -q ${files} \
        | awk -v OFS='\t' \
            -v val='${group_key}' \
            '
                NR==1 {print \$0, "group_id"}
                NR>1 {print \$0, val}
            '
        | sort-bed - >> ${name}
    """
}



process pack_data {
    conda params.conda
    publishDir params.outdir
    scratch true
    label "high_mem"

    input:
        path aggregated_variants
        path non_aggregated_variants
        val aggregation_key
    
    output:
        path sorted_aggregated, emit: aggregated
        tuple path(sorted_non_aggregated), path("${sorted_non_aggregated}.tbi"), emit: non_aggregated
        tuple path(name), path(parquet_non_aggregated), emit: stats
    
    script:
    sorted_aggregated = "aggregated.${aggregation_key}.bed"
    sorted_non_aggregated = "non_aggregated.${aggregation_key}.bed.gz"
    parquet_non_aggregated = "non_aggregated.${aggregation_key}.parquet"
    name = "cav_stats.${aggregation_key}.tsv"
    """
    head -1 ${aggregated_variants} > ${sorted_aggregated}
    sort-bed ${aggregated_variants} >> ${sorted_aggregated}

    head -1 ${non_aggregated_variants} > tmp.txt
    sort-bed ${non_aggregated_variants} >> tmp.txt
    bgzip -c tmp.txt > ${sorted_non_aggregated}
    tabix ${sorted_non_aggregated}

    python3 -c "import pandas as pd; pd.read_table('${sorted_non_aggregated}').to_parquet('${parquet_non_aggregated}')"

    python3 $moduleDir/bin/collect_stats.py \
        ${sorted_aggregated} \
        ${name} \
        --fdr ${params.fdr_tr}
    """

}

process get_correspondence_map {

    conda params.conda

    input:
        path samples_file

    output:
        path name

    script:
    name = "${params.aggregation_key}.correspondence_map.tsv"
    """
    python3 $moduleDir/bin/get_value_for_agg_key.py \
        ${samples_file} \
        ${params.aggregation_key} \
        ${name}
    """
}

workflow packData {
    take:
        merged_files
        aggregation_key

    main:
        aggregated_merged = merged_files
            | aggregate_pvals
            | map(it -> it[1])
            | collectFile(
                skip: 1,
                keepHeader: true,
                name: "aggregated.not_sorted.${aggregation_key}.bed",
            )
        non_aggregated_merged = merged_files
            | map(it -> it[1])
            | collectFile(
                skip: 1,
                keepHeader: true,
                name: "non_aggregated.not_sorted.${aggregation_key}.bed",
            )
        
        packed = pack_data(aggregated_merged, non_aggregated_merged, aggregation_key)
    emit:
        packed.aggregated
        packed.non_aggregated
}

workflow aggregation {
    take:
        sample_split_pvals
    main:
        params.aggregation_key = params.aggregation_key ?: "all"
        sample_cl_correspondence = Channel.fromPath(params.samples_file)
            | get_correspondence_map
            | splitCsv(header: true, sep: '\t')
            | map(row -> tuple(row.sample_id, row.value))

        merged_pvals = sample_split_pvals // sample_id, file
            | join(sample_cl_correspondence) // sample_id, file, value
            | filter(it -> !it[2].isEmpty())
            | map(it -> tuple(it[2], it[1])) // value, file
            | set_key_for_group_tuple 
            | groupTuple()
            | merge_files

        packed = packData(
            merged_pvals,
            params.aggregation_key
        )
    emit:
        packed[0]
        packed[1]
}


///////////////////////////////////////////////////////////////
/////////////////////// Entry workflows ///////////////////////
///////////////////////////////////////////////////////////////
workflow {
    Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))
        | aggregation
}

process filter_bad1 {
    conda params.conda
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(non_aggregated)
    
    output:
        tuple val(sample_id), path(name)
    script:
    name = "${sample_id}_BAD1.nonaggregated.bed"
    """
    awk -v OFS='\t' -F'\t' '
        NR == 1 { \
            for (i=1; i<=NF; i++) { \
                if (\$i == "BAD") { \
                    badCol = i; \
                    print \$0; \
                    break; \
                } \
            } \
        } \
        NR > 1 && \$(badCol) == 1.0 { \
            print \$0; \
        }' ${non_aggregated} > ${name}
    """
}

workflow bad1Aggregation {
    // This workflow is used to do overall aggregation for variants in bad1 regions
    // It is agnostic to params.aggregation_key
    pvals = Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))
        | filter_bad1
        | map(it -> it[1])
        | collect(sort: true)
        | map(it -> tuple('all.bad1', it)) 
        | merge_files
    packData(pvals, 'all.bad1')
}

// DEFUNC
workflow tmp {
    Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))
        | filter_bad1
        | aggregation
}