#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process extract_indiv_vcfs {
    tag "${indiv_id}"
    conda params.conda
    publishDir "${params.outdir}/filtered_vcfs"

    input:
	    tuple val(indiv_id), val(agg_numbers)
    output:
        tuple val(indiv_id), path(name), path("${name}.csi")
    script:
    name = "${indiv_id}.vcf.gz"
    """
    bcftools view --output-type z -s ${agg_numbers} ${params.vcf_file} > ${name}
    bcftools index ${name}
    """
}

process extract_ag_id_vcf {
    tag "${ag_id}"
    conda params.conda
    publishDir "${params.outdir}/filtered_vcfs"

    input:
	    val ag_id
    output:
        tuple val(ag_id), path(name), path("${name}.csi")
    script:
    name = "${ag_id}.bed"
    """
    bcftools query -s ${ag_id} -i'GT="alt"' \
      -f'%CHROM\\t%POS0\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%AD{0}\\t%AD{1}\\t%GT]\\t%INFO/TOPMED\n' \
        ${params.vcf_file} > ${name}
    """
}

// BABACHI filter
process filter_indiv_vcfs {
    tag "${indiv_id}"
    publishDir "${params.outdir}/filtered_bed"
    conda params.conda

    input:
	    tuple val(indiv_id), path(indiv_vcf)
    output:
        tuple val(indiv_id), path(name)
    script:
    name = "${indiv_id}.snps.bed"
    """
    babachi filter ${indiv_vcf} -O ${name} -a ${params.allele_tr}
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

workflow test {
    ag_ids = Channel.of('AG3864', 'AG3831', 'AG3813', 'AG3958', 'AG32073')
    extract_ag_id_vcf(ag_ids)
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
        extractAggNumbers(sample_ag_merge) | filter_indiv_vcfs
    emit:
        filter_indiv_vcfs.out
}


workflow {
    extractAndFilter()
}