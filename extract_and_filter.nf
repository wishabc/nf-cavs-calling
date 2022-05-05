#!/usr/bin/env nextflow

process extract_indiv_vcfs {
    tag "Extracting ${indiv_id}"
    publishDir params.outdir + '/indiv_vcfs'

    input:
	    tuple val(indiv_id), val(agg_numbers)
    output:
        tuple val(indiv_id), path("${indiv_id}.vcf")
    script:
    """
    bcftools view -s ${agg_numbers} ${params.vcfFile} > ${indiv_id}.vcf
    """
}
workflow extract_agg_numbers {
    take:
        tuple id_agg_numbers
    main:
        extract_indiv_vcfs(id_agg_numbers)
    emit:
        extract_indiv_vcfs.out
}


process filter_indiv_vcfs {
    tag "Filtering ${indiv_id}"
    publishDir params.outdir + 'filtered_indiv_vcfs'

    input:
	    tuple val(indiv_id), path(indiv_vcf)
    output:
        tuple val(indiv_id), path(name)
    script:
    name = get_filtered_file_by_indiv_id(indiv_id)
    """
    babachi filter ${indiv_vcf} -O ${name}
    rm ${indiv_vcf}
    """
}

def get_filtered_file_by_indiv_id(ind) {
    "${ind}.snps.bed"
}


workflow extractAndFilter {
    main:
        sample_ag_merge = Channel
                .fromPath(params.samplesFile)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id, row.ag_number) }
                .groupTuple(by:0)
                .map{ it -> tuple(it[0], it[1].join(",")) }
        extract_agg_numbers(sample_ag_merge) | filter_indiv_vcfs 
    emit:
        filter_indiv_vcfs.out
}

workflow extractAggNumbers {
    main:
        ag_merge = Channel
                .fromPath(params.samplesFile)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.ag_number, row.ag_number)}
        extract_agg_numbers(ag_merge)
    out:
        extract_agg_numbers.out   
}

workflow {
    extractAndFilter()
}

