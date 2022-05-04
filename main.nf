#!/usr/bin/env nextflow

params.mode = true
params.samples_file = 'metadata+ag_number.tsv'
params.outdir = ''
params.states = "1,1.5,2,2.5,3,4,5,6"
params.prior = "uniform"

// DO NOT EDIT BELOW

Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple(row.indiv_id, row.ag_number) }
	.groupTuple(by:0)
	.map{ it -> tuple(it[0], it[1].join(",")) }
	.set{ sample_ag_merge }


process extract_indiv_vcfs {
    tag "${indiv_id}"
    publishDir params.outdir + '/filtered_indiv_vcfs'

    input:
	    tuple val(indiv_id), val(agg_numbers) from sample_ag_merge
    output:
        tuple val(indiv_id), path("${indiv_id}.snps.bed") into indiv_snps_file
    script:
    """
    echo 'STARTING $indiv_id'
    bcftools view -s ${agg_numbers} allele_counts.fixed.vcf.gz > ${indiv_id}.vcf
    babachi filter ${indiv_id}.vcf -O ${indiv_id}.snps.bed
    """
}

