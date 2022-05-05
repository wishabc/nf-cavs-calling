#!/usr/bin/env nextflow

process extract_indiv_vcfs {
    tag "Extracting ${indiv_id}"
    publishDir params.outdir + '/indiv_vcfs'

    input:
	    tuple val(indiv_id), val(agg_numbers)
    output:
        tuple val(indiv_id), path("${indiv_id}.vcf.gz"), path("${indiv_id}.vcf.gz.csi")
    script:
    """
    bcftools view -s ${agg_numbers} ${params.vcfFile} > ${indiv_id}.vcf
    bgzip ${indiv_id}.vcf
    tabix ${indiv_id}.vcf.gz
    """
}

// BABACHI filter
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
                .fromPath(params.samplesFile)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id, row.ag_number) }
                .groupTuple(by:0)
                .map{ it -> tuple(it[0], it[1].join(",")) }
        extractAggNumbers(sample_ag_merge) | filter_indiv_vcfs 
    emit:
        filter_indiv_vcfs.out
}


workflow {
    extractAndFilter()
}

// Extract each sample in separate file
workflow extractAllSamples {
    main:
        ag_merge = Channel
                .fromPath(params.samplesFile)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id + '@' + row.ag_number, row.ag_number)}
        extractAggNumbers(ag_merge)
    emit:
        extractAggNumbers.out
}
