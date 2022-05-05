#!/usr/bin/env nextflow

include {get_filtered_file_by_indiv_id} from "./extract_and_filter"

process apply_babachi {
	cpus 2
    tag "BABACHI ${indiv_id}"
    publishDir params.outdir + '/badmaps'

	input:
		tuple val(indiv_id), path(snps_file)
	output:
		tuple val(indiv_id), path("${indiv_id}.bad.bed")

	script:
	"""
    babachi ${snps_file} -O ${indiv_id}.bad.bed --visualize -z -e png -j ${task.cpus} -p ${params.prior} -s ${params.states}
	"""
}


process intersect_with_snps {

    publishDir params.outdir + '/snp_annotation'
	input:
		tuple val(indiv_id), path(snps_file), path(badmap_file)
    
    output:
        tuple val(indiv_id), path("${indiv_id}.intersect.bed")

	script:
	"""
	bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb | awk \
    '{for (i=1;i<=7;i+=1) print \$i, print \$11}' > ${indiv_id}.intersect.bed
	"""
}


def get_filtered_vcf_path(filtered_vcf_path, indiv_id) {
    file = get_filtered_file_by_indiv_id(indiv_id)
    if (filtered_vcf_path != '') {
        return "${filtered_vcf_path}/${file}"
    }
    else {
        return file
    }
}

workflow estimate_bad {
    main:
        extracted_vcfs = Channel.fromPath(params.samplesFile)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.indiv_id,
                get_filtered_vcf_path(params.filteredVcfs, row.indiv_id)))
            .distinct()
        extracted_vcfs.last().view()
        badmaps_map = apply_babachi(extracted_vcfs)
        badmaps_and_snps = extracted_vcfs.join(
            badmaps_map
        )
        intersect_with_snps(badmaps_and_snps)
    emit:
        intersect_with_snps.out
}

workflow {
    estimate_bad()
}
