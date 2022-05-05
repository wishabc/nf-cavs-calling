#!/usr/bin/env nextflow

include {get_filtered_file_by_indiv_id} from "./extract_and_filter"

process apply_babachi {
	cpus 2
    publishDir params.outdir + '/badmaps'

	input:
		tuple val(indiv_id), path(snps_file)
	output:
		tuple val(indiv_id), path(snps_file), path("${indiv_id}.bad.bed")

	script:
	"""
    babachi ${snps_file} -O ${indiv_id}.bad.bed --visualize -z -e png \
	 -j ${task.cpus} -p ${params.prior} -s ${params.states}
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

workflow estimate_bad {
    main:
        extracted_vcfs = Channel.fromPath(params.samplesFile)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.indiv_id,
                params.filteredVcfs + '/' + get_filtered_file_by_indiv_id(row.indiv_id)) }

        apply_babachi(extracted_vcfs) | intersect_with_snps
    emit:
        intersect_with_snps.out
}

workflow {
    estimate_bad()
}
