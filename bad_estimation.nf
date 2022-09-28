#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process apply_babachi {
	cpus 2
    tag "${indiv_id}"
    publishDir "${params.outdir}/${outpath}badmaps", pattern: "${badmap_file}"
    publishDir "${params.outdir}/${outpath}intersect", pattern: "${name}"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file)
        val outpath
	output:
		tuple val(indiv_id), path(name)

	script:
    badmap_file = "${indiv_id}.bad.bed"
    name = "${indiv_id}.${outpath}intersect.bed"
	"""
    babachi ${snps_file} -O ${badmap_file} -j ${task.cpus} -p ${params.prior} -s ${params.states} -a ${params.allele_tr} --geometric-prior ${params.geometric_prior}
    head -1 ${badmap_file} | xargs -I % echo "#chr\tstart\tend\tID\tref\talt\tref_counts\talt_counts\tsample_id\t%" > ${name}
    if [[ \$(wc -l <${snps_file}) -ge 2 ]]; then
	    bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb >> ${name}
    fi
	"""
}

workflow estimateBad {
    take:
        extracted_vcfs
        outpath
    main:
        out = apply_babachi(extracted_vcfs
            .filter { it[1].countLines() > 1 }, outpath)
            .filter { it[1].countLines() > 1 }
    emit:
        out
}

workflow estimateBadByIndiv {
    take:
        prefix
    main:
        filtered_vcfs = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.indiv_id, file(row.snps_file)))

        badmaps_map = estimateBad(filtered_vcfs, prefix) 
    emit:
        filtered_vcfs
        out
}

workflow {
    estimateBadByIndiv('')
}