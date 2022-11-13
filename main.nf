#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; addImputedCavs; aggregate_pvals } from "./pval_calc"
include { motifEnrichment } from "./motif_enrichment"


params.conda = "$moduleDir/environment.yml"


process sort_and_gzip {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path inp

    output:
        tuple path(name), path("${name}.tbi")

    script:
    name = "${inp.simpleName}.sorted.bed.gz"
    """
    cat ${inp} | grep -v '^#' | sort-bed - | bgzip -c > ${name}
    tabix ${name}
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


workflow test {
    binom_p = Channel.fromPath('/net/seq/data2/projects/sabramov/ENCODE4/cav-calling/babachi_1.5_common_final/output/final.pval_files_binom/*.bed')
        .map(it -> tuple(it.simpleName, file(it)))
    sample_split_pvals = split_into_samples(binom_p)
        .flatten()
        .map(it -> tuple(it.simpleName, it))
    
    if (params.aggregation_key && params.aggregation_key != 'all') {
        sample_cl_correspondence = Channel.fromPath(params.samples_file)
                .splitCsv(header:true, sep:'\t')
                .map(row -> tuple(row.ag_id, row[params.aggregation_key]))
        pvals = sample_split_pvals
        .join(sample_cl_correspondence)
        .collectFile(keepHeader: true, skip: 1) { item -> [ "${item[2]}.bed", item[1].text + '\n' ]}
        .map(it -> tuple(it.simpleName, it))
    } else {
        pvals = sample_split_pvals.collectFile()
        .map(it -> tuple('all', it))
    }

    aggregate_pvals(pvals, 'binom', 'final.')  // | map(it -> it[1]) | motifEnrichment
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
    binom_p = calcPvalBinom(imputed_cavs, iter2_prefix)
    sample_split_pvals = split_into_samples(binom_p)
        .flatten()
        .map(it -> tuple(it.simpleName, it))
    
    if (params.aggregation_key && params.aggregation_key != 'all') {
        sample_cl_correspondence = Channel.fromPath(params.samples_file)
                .splitCsv(header:true, sep:'\t')
                .map(row -> tuple(row.ag_id, row[params.aggregation_key]))
        pvals = sample_split_pvals
        .join(sample_cl_correspondence)
        .collectFile(keepHeader: true, skip: 1) { item -> [ "${item[2]}.bed", item[1].text + '\n' ]}
        .map(it -> tuple(it.simpleName, it))
    } else {
        pvals = sample_split_pvals.collectFile(
                name: "all_variants.bed"
        )
        .map(it -> tuple('all', it))
    }

        // .concat(all_pval_file)
    aggregate_pvals(pvals, 'binom', 'final.')
}