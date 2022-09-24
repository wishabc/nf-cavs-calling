#!/usr/bin/env nextflow
include { get_file_by_indiv_id; get_id_by_sample } from "./helpers"

params.conda = "$moduleDir/environment.yml"

process apply_babachi {
	cpus 2
    tag "BABACHI ${indiv_id}"
    publishDir "${params.outdir}/${outpath}badmaps"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file)
        val outpath
	output:
		tuple val(indiv_id), path(name)

	script:
    name = get_file_by_indiv_id(indiv_id, "badmap")
	"""
    babachi ${snps_file} -O ${name} -j ${task.cpus} -p ${params.prior} -s ${params.states} -a ${params.alleleTr} --geometric-prior ${params.geometricPrior}
	"""
}


process intersect_with_snps {
    tag "Annotating SNPs ${indiv_id}"
    publishDir "${params.outdir}/${outpath}intersect"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file), path(badmap_file)
        val outpath
    output:
        tuple val(indiv_id), path(name)

	script:
    name = get_file_by_indiv_id(indiv_id, "${outpath}intersect")
	"""
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
        out = apply_babachi(extracted_vcfs.filter { it[1].countLines() > 1 }, outpath).filter { it[1].countLines() > 1 }
    emit:
        out
}

workflow intersectWithBadmap {
    take:
        badmaps_and_snps
        outpath
    main:
        out = intersect_with_snps(badmaps_and_snps, outpath).filter { it[1].countLines() > 1 }
    emit:
        out
}

workflow estimateBadByIndiv {
    main:
        filtered_vcfs = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.indiv_id, file(row.snps_file)))

        badmaps_map = estimateBad(filtered_vcfs, '') 
        badmaps_and_snps = filtered_vcfs.join(
            badmaps_map
        )
        out = intersectWithBadmap(badmaps_and_snps, '')

    emit:
        filtered_vcfs
        out
}

workflow {
    estimateBadByIndiv()
}