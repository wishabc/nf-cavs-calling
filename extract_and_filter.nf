#!/usr/bin/env nextflow
include { get_file_by_indiv_id } from "./helpers"

params.conda = "$moduleDir/environment.yml"

process extract_indiv_vcfs {
    tag "${indiv_id}"
    conda params.conda

    input:
	    tuple val(indiv_id), val(agg_numbers)
    output:
        tuple val(indiv_id), path(name), path("${name}.csi")
    script:
    name = get_file_by_indiv_id(indiv_id, "vcf")
    """
    bcftools view --output-type z -s ${agg_numbers} ${params.vcfFile} > ${name}
    bcftools index ${name}
    """
}

// BABACHI filter
process filter_indiv_vcfs {
    tag "${indiv_id}"
    publishDir "${params.outdir}/filtered_vcfs"
    conda params.conda

    input:
	    tuple val(indiv_id), path(indiv_vcf)
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_file_by_indiv_id(indiv_id, "filter")
    """
    babachi filter ${indiv_vcf} -O ${name} -a ${param.alleleTr}
    """
}

process extend_metadata {
    publishDir params.outdir
    conda params.conda

	input:
		tuple val(indiv_id), path(bed_files)

	output:
		path name

	script:
	name = 'metadata+filtered_vcfs.txt'
	column = ['filtered_vcf', *bed_files].join('\n')
	"""
	echo "${column}" > columns.txt
	paste ${params.samples_file} columns.txt > ${name}
	"""
}

// Extract samples by map file
workflow extractAggNumbers {
    take: id_agg_numbers
    main:
        extract_indiv_vcfs(id_agg_numbers)
    emit:
        extract_indiv_vcfs.out.map(
            it -> tuple(it[0], it[1])
        )
}

// Extract samples grouped by INDIV_ID
workflow extractAndFilter {
    main:
        sample_ag_merge = Channel
                .fromPath(params.samples_file)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id, row.ag_number) }
                .groupTuple(by:0)
                .map{ it -> tuple(it[0], it[1].join(",")) }
        extractAggNumbers(sample_ag_merge) | filter_indiv_vcfs | extend_metadata
    emit:
        filter_indiv_vcfs.out
}


workflow {
    extractAndFilter()
}

// // Defunc
// Extract each sample in separate file
// workflow extractAllSamples {
//     main:
//         ag_merge = Channel
//                 .fromPath(params.samples_file)
//                 .splitCsv(header:true, sep:'\t')
//                 .map(row -> tuple(get_id_by_sample(row.indiv_id, row.ag_number), row.ag_number))
//         extractAggNumbers(ag_merge)
//     emit:
//         extractAggNumbers.out
// }
// 
// // Filter each sample with BABACHI
// workflow filterAllSamples {
//     main:
//         ag_merge = Channel
//                 .fromPath(params.samples_file)
//                 .splitCsv(header:true, sep:'\t')
//                 .map(row -> get_id_by_sample(row.indiv_id, row.ag_number))
//                 .map(it -> tuple(it,
//                  raw_vcfs_dir + get_filtered_file_by_indiv_id(it, 'vcf')))
//         filter_indiv_vcfs(ag_merge)
//     emit:
//         filter_indiv_vcfs.out

// }

// workflow extractAndFilterAllSamples {
//     extractAllSamples | filterAllSamples
// }