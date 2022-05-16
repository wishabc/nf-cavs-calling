#!/usr/bin/env nextflow

include { get_file_by_indiv_id; get_id_by_sample } from "./helpers"


def get_filtered_vcf_path(filtered_vcf_path, indiv_id) {
    file = get_file_by_indiv_id(indiv_id, "filter")
    if (filtered_vcf_path != '') {
        return "${filtered_vcf_path}/${file}"
    }
    else {
        return file
    }
}

process apply_babachi {
	cpus 2
    tag "BABACHI ${indiv_id}"
    publishDir params.outdir + '/badmaps'

	input:
		tuple val(indiv_id), path(snps_file)
	output:
		tuple val(indiv_id), path(name)

	script:
    name = get_file_by_indiv_id(indiv_id, "badmap")
	"""
    babachi ${snps_file} -O ${name} --visualize -z -e png -j ${task.cpus} -p ${params.prior} -s ${params.states}
	"""
}


process intersect_with_snps {
    tag "Annotating SNPs ${indiv_id}"
    publishDir params.outdir + '/snp_annotation'
	input:
		tuple val(indiv_id), path(snps_file), path(badmap_file)
    
    output:
        tuple val(indiv_id), path(name)

	script:
    name = get_file_by_indiv_id(indiv_id, "intersect")
	"""
    head -1 ${badmap_file} | xargs -I % echo "#chr\tstart\tend\tID\tref\talt\tref_counts\talt_counts\t%" > ${name}
	bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb >> ${name}
	"""
}


workflow estimateBadAndIntersect {
    take:
        extracted_vcfs
    main:
        badmaps_map = apply_babachi(extracted_vcfs)
        badmaps_and_snps = extracted_vcfs.join(
            badmaps_map
        )
        intersect_with_snps(badmaps_and_snps)
    emit:
        intersect_with_snps.out
}

workflow estimateBadByIndiv {
    main:
        filtered_vcfs = Channel.fromPath(params.samplesFile)
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.indiv_id,
            get_filtered_vcf_path(params.filteredVcfs, row.indiv_id)))
        .distinct()
        estimateBadAndIntersect(filtered_vcfs)
    emit:
        estimateBadAndIntersect.out
}

workflow {
    if (params.filteredVcfs == "") {
        params.filteredVcfs = "${params.outdir}/filtered_vcfs"
    }
    estimateBadByIndiv()
}

workflow estimateBadBySample {
    filtered_vcfs = Channel.fromPath(params.samplesFile)
    .splitCsv(header:true, sep:'\t')
    .map(row -> get_id_by_sample(row.indiv_id, row.ag_number))
    .map(it -> tuple(it,
        get_filtered_vcf_path(params.filteredVcfs, it)))
    estimateBadAndIntersect(filtered_vcfs)
}