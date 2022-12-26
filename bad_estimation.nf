#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process apply_babachi {
	cpus 2
    scratch true
    tag "${indiv_id}"
    publishDir "${params.outdir}/${outpath}badmaps", pattern: "${badmap_file}"
    publishDir "${params.outdir}/${outpath}intersect", pattern: "${name}"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file)
        val outpath
	output:
		tuple val(indiv_id), path(name), emit: intersect
        tuple val(indiv_id), path(badmap_file), emit: babachi

	script:
    badmap_file = "${indiv_id}.bad.bed"
    name = "${indiv_id}.${outpath}intersect.bed"
    prior_params = params.prior == 'geometric' ? "--geometric-prior ${params.geometric_prior}" : ""
	"""
    cat ${snps_file} | \$10 >= 0.05 {print;} > snps.common.bed
    if [[ `wc -l < snps.common.bed` -le 100 ]]; then
	    touch ${name}
        touch ${badmap_file}
        exit 0
    fi
    babachi snps.common.bed -O ${badmap_file} -j ${task.cpus} \
        -p ${params.prior} ${prior_params} \
        -s ${params.states} -a ${params.allele_tr}


    head -1 ${badmap_file} | xargs -I % echo "#chr\tstart\tend\tID\tref\talt\tref_counts\talt_counts\tsample_id\tMAF\tFMR\t%" > ${name}
    
    # Avoid intersecting with empty file
    if [[ `wc -l <${snps_file}` -ge 2 ]]; then
	    bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb >> ${name}
    fi
	"""
}

workflow estimateBad {
    take:
        extracted_snps
        outpath
    main:
        non_empty_snps = extracted_snps.filter { it[1].countLines() >= 100 }
        out = apply_babachi(non_empty_snps, outpath).intersect
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
            .map(row -> tuple(row.indiv_id, file(row.snps_file, )))
            .filter { it[1].exists() }
            .unique()

        out = estimateBad(filtered_vcfs, prefix) 
    emit:
        filtered_vcfs
        out
}

workflow {
    estimateBadByIndiv('')
}