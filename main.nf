#!/usr/bin/env nextflow
include { callCavsFromVcfsBinom; calcPvalBinom; addExcludedCavs; aggregate_pvals } from "./pval_calc"

params.conda = "$moduleDir/environment.yml"


def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  | transpose()
}


process apply_babachi {
	cpus 3
    scratch true
    tag "${indiv_id}"
    publishDir "${params.outdir}/${prefix}.badmaps", pattern: "${badmap_file}"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file)
        val prefix
	output:
		tuple val(indiv_id), path(name), emit: intersect
        tuple val(indiv_id), path(badmap_file), emit: babachi

	script:
    badmap_file = "${indiv_id}.bad.bed"
    name = "${indiv_id}.${prefix}.intersect.bed"
    prior_params = params.prior == 'geometric' ? "--geometric-prior ${params.geometric_prior}" : ""
	"""
    head -n 1 ${snps_file} > header.txt
    tail -n +2 ${snps_file} | awk '((\$10 >= 0.05) && (\$10 != "None")) { print; }' > snps.common.bed
    if [[ `wc -l < snps.common.bed` -le ${params.min_snps_count} ]]; then
	    touch ${name}
        touch ${badmap_file}
        exit 0
    fi

    babachi snps.common.bed \
        -O ${badmap_file} \
        -j ${task.cpus} \
        -p ${params.prior} \
        ${prior_params} \
        -s ${params.states} \
        -a ${params.allele_tr}


    head -1 ${badmap_file} | xargs -I % echo "`cat header.txt`\t%" > ${name}
    
    # Avoid intersecting with empty file
    if [[ `wc -l <${snps_file}` -ge 2 ]]; then
	    bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb >> ${name}
    fi
	"""
}


process merge_files {
    conda params.conda
    publishDir "${params.outdir}/nonaggregated.${params.aggregation_key}", pattern: "${name}"
    scratch true
    tag "${group_key}"

    input:
        tuple val(group_key), path(files)

    output:
        tuple val(group_key), path(name)

    script:
    name = "${group_key}.sorted.bed"
    """
    python3 $moduleDir/bin/merge_files.py f.txt ${group_key} ${files}
    head -n1 f.txt > "${name}"
    sort-bed f.txt >> "${name}"
    """
}


process split_into_samples {
    tag "${indiv_id}"
    conda params.conda

    input:
        tuple val(indiv_id), path(pval_file)
    
    output:
        path "*${suffix}"
    
    script:
    suffix = "sample_split.bed"
    """
    python3 $moduleDir/bin/split_into_samples.py ${pval_file} ${suffix}
    """
}


process annotate_variants {
    conda params.conda
    tag "${sample_id}"
    publishDir "${params.outdir}/by_sample"
    scratch true

    input:
        tuple val(sample_id), path(pval_file), val(hotspots_file), val(footprint_file)

    output:
        tuple val(sample_id), path(name)

    script:
    name = "${sample_id}.nonaggregated.bed"
    """
    sort-bed ${pval_file} > pval_f.bed

    # Add footprints
    if [[ "${footprint_file}" == "" ]]; then
        cat pval_f.bed | awk '{ print "-" }' > footprints.txt
    else
        bedmap --header \
            --indicator pval_f.bed \
            ${footprint_file} > footprints.txt
    fi

    if [[ "${hotspots_file}" == "" ]]; then
        cat pval_f.bed | awk '{ print "-" }' > hotspots.txt
    else
        bedmap --header \
            --indicator pval_f.bed \
            ${hotspots_file} > hotspots.txt
    fi

    echo -e "`head -1 ${pval_file}`\tfootprints\thotspots" > ${name}
    paste pval_f.bed footprints.txt hotspots.txt >> ${name}
    """
}

workflow estimateBAD {
    take:
        snv_files
        prefix
    main:
        non_empty = snv_files
            | filter { it[1].countLines() >= params.min_snps_count }
        out = apply_babachi(non_empty, prefix).intersect
            | filter { it[1].countLines() > 1 } // At least one line except header
    emit:
        out
}

workflow estimateBADByIndiv {
    take:
        prefix
    main:
        filtered_vcfs = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.indiv_id, file(row.snps_file)))
            .unique()
        out = estimateBAD(filtered_vcfs, prefix) 
    emit:
        filtered_vcfs
        out
}

workflow aggregation {
    take:
        sample_split_pvals
    main:
        params.aggregation_key = params.aggregation_key ?: "all"
        if (params.aggregation_key != 'all') {
            sample_cl_correspondence = Channel.fromPath(params.samples_file)
                    | splitCsv(header:true, sep:'\t')
                    | map(row -> tuple(row.ag_id, row[params.aggregation_key]))
    
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

        merged = merge_files(pvals)
        out = aggregate_pvals(merged, iter2_prefix)
        
        if (agg_key != "all") {
            out.collectFile(
                storeDir: params.outdir,
                name: "aggregated.${agg_key}.bed",
            )
            non_aggregated_merge = merged.collectFile(
                storeDir: params.outdir,
                name: "non_aggregated.${agg_key}.bed",
            )
        } else {
            non_aggregated_merge = merged
        }
    emit:
        out
        non_aggregated_merge
}

workflow annotateWithFootprints {
    take:
        pval_files
    main:
        annotations = Channel.fromPath(params.samples_file)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.ag_id,
                    row?.hotspots_file ? file(row.hotspots_file) : null, 
                    row?.footprints_file ? file(row.footprints_file) : null)
                )
        out = pval_files
            | join(annotations)
            | annotate_variants
    emit:
        out
}


workflow {
    // Estimate BAD and call 1-st round CAVs
    iter1_prefix = 'iter1'

    bads = Channel.of(params.states.tokenize(','))
    filtered_vcfs_and_intersect = estimateBADByIndiv(iter1_prefix)
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files, iter1_prefix)

    iter2_prefix = 'final'
    // Reestimate BAD, and add excluded SNVs
    all_snps = estimateBAD(no_cavs_snps, iter2_prefix)
        | join(intersect_files)
        | addExcludedCavs

    // Annotate with footprints and hotspots + aggregate by provided aggregation key
    agg_files = calcPvalBinom(all_snps, iter2_prefix)
        | split_into_samples
        | flatten()
        | map(it -> tuple(it.simpleName, it))
        | annotateWithFootprints
        | aggregation
}   


// Aggregation only workflow
workflow aggregatePvals {
    Channel.fromPath("$launchDir/${params.outdir}/by_sample/*.bed")
        | map(it -> tuple(it.simpleName, it))
        | aggregation
}
