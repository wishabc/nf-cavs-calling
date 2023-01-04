#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; addImputedCavs; aggregate_pvals } from "./pval_calc"

def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}


params.conda = "$moduleDir/environment.yml"


process merge_files {
    conda params.conda
    publishDir "${params.outdir}/pvals_nonaggregated.${params.aggregation_key}", pattern: "${name}"
    scratch true
    tag "${group_key}"

    input:
        tuple val(group_key), path(files)

    output:
        tuple val(group_key), path(name)

    script:
    name = "${group_key}.sorted.bed"
    """
    python3 $moduleDir/bin/merge_files.py f.txt ${files}
    head -n1 f.txt > ${name}
    sort-bed f.txt >> ${name}
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
    publishDir "${params.outdir}/annotations"
    scratch true

    input:
        tuple val(sample_id), path(pval_file), path(hotspots_file), val(footprint_file)

    output:
        tuple val(sample_id), path(name)

    script:
    name = "${sample_id}.fp_annotation.bed"
    footprint_f = footprint_file ? footprint_file : "empty.bed"
    """
    if [[ "${footprint_f}" == "empty.bed" ]]; then
        touch ${footprint_f}
    fi
    sort-bed ${pval_file} > pval_f.bed

    bedmap --header \
        --indicator pval_f.bed \
        ${footprint_f} >> footprints.txt
    
    bedmap --header \
        --indicator pval_f.bed \
        ${hotspots_file} >> hotspots.txt

    echo -e "`head -1 ${pval_file}`\tfootprints\thotspots" > ${name}
    paste pval_f.bed footprints.txt hotspots.txt >> ${name}
    """
}

workflow aggregation {
    take:
        sample_split_pvals
    main:
        agg_key = params.aggregation_key ? params.aggregation_key : "all"
        if (agg_key != 'all') {
            sample_cl_correspondence = Channel.fromPath(params.samples_file)
                    .splitCsv(header:true, sep:'\t')
                    .map(row -> tuple(row.ag_id, row[params.aggregation_key]))
            pvals = sample_split_pvals
                | join(sample_cl_correspondence)
                | filter(it -> !it[2].isEmpty())
                | map(it -> tuple(it[2], it[1]))
                | set_key_for_group_tuple 
                | groupTuple() 
        } else {
            pvals = sample_split_pvals
                | map(it -> it[1])
                | collect()
                | map(it -> tuple('all', it)) 
                

        }
        out = aggregate_pvals(merge_files(pvals), "binom.${agg_key}", 'final.')
    emit:
        out
}

workflow annotateWithFootprints {
    take:
        pval_files
    main:
        annotations = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.ag_id,
                                file(row.hotspots_file), 
                                row?.footprint_path ? file(row.footprint_path) : null)
                )
        data = pval_files.join(annotations)
        annotations = annotate_variants(data)
    emit:
        annotations
}


workflow {
    // Estimate BAD and call 1-st round CAVs
    iter1_prefix = 'iter1.'

    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv(iter1_prefix)
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files, iter1_prefix)

    iter2_prefix = 'final.'
    // Reestimate BAD, and add excluded SNVs
    iter2_intersections = estimateBad(no_cavs_snps, iter2_prefix)
    imputed_cavs = addImputedCavs(iter2_intersections.join(intersect_files))
    binom_p = calcPvalBinom(imputed_cavs, iter2_prefix)[0]
    sample_split_pvals = split_into_samples(binom_p).flatten()
        .map(it -> tuple(it.simpleName, it))

    // Annotate with footprints and hotspots + aggregate by provided aggregation key
    ann_pvals = annotateWithFootprints(sample_split_pvals)
    aggregation(ann_pvals)
}


// Only aggregation workflow
params.sample_pvals_dir = "$launchDir/${params.outdir}/annotations"
workflow aggregatePvals {
    sample_pvals = Channel.fromPath("${params.sample_pvals_dir}/*.bed")
        .map(it -> tuple(file(it).simpleName, file(it)))
    aggregation(sample_pvals)
}


workflow withFootprints {
    sample_pvals = Channel.fromPath("${params.sample_pvals_dir}/*.bed")
        .map(it -> tuple(file(it).simpleName, file(it)))
    
    annotateWithFootprints(sample_pvals)
}