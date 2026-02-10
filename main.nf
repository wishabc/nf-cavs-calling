#!/usr/bin/env nextflow
include { callCavsFirstRound; calcPvalBinom; add_cavs } from "./pval_calc"
include { aggregation } from "./aggregation"

params.conda = "$moduleDir/environment.yml"


def check_var(var, prefix) {
    if (var) {
        if (file(var).exists()) {
            return file(var)
        }
    } 
    return file("${prefix}.empty") 
}

process collect_files {
    conda params.conda
    tag "${badmap_id}"

    input:
        tuple val(badmap_id), path(babachi_files)

    output:
        tuple val(badmap_id), path(name)
    
    script:
    name = "${badmap_id}.bed"
    """
    head -n 1 ${babachi_files[0]} > ${name}
    tail -n +2 -q ${babachi_files} \
        | awk '\$7+\$8 >= ${params.initial_coverage_filter} {print;}' \
        | sort-bed - >> ${name}
    """
}



process apply_babachi {
	cpus 3
    scratch true
    tag "${indiv_id}"
    publishDir "${params.outdir}/${prefix}.badmaps", pattern: "${badmap_file}"
    conda params.conda

	input:
		tuple val(indiv_id), path(snps_file)
        val prefix
	output:
		tuple val(indiv_id), path(name), emit: intersect
        tuple val(indiv_id), path(badmap_file), emit: babachi

	script:
    badmap_file = "${indiv_id}.bad.bed"
    name = "${indiv_id}.${prefix}.intersect.bed"
    prior_params = params.prior == 'geometric' ? "--geometric-prior ${params.geometric_prior}" : ""
	"""
    head -n 1 ${snps_file} > header.txt
    tail -n +2 ${snps_file} | awk '( \
        (\$10 >= ${params.babachi_maf_tr}) \
            && (\$10 != "None") \
            && (\$11 >= ${params.babachi_maf_tr}) \
            && (\$11 != "None") \
    ) { print; }' > snps.common.bed
    if [[ `wc -l < snps.common.bed` -le ${params.min_snps_count} ]]; then
	    touch ${name}
        touch ${badmap_file}
        exit 0
    fi

    # -j ${task.cpus} 
    babachi snps.common.bed \
        -O ${badmap_file} \
        -p ${params.prior} \
        ${prior_params} \
        -s ${params.states} \
        -a ${params.babachi_allele_tr}


    head -1 ${badmap_file} | xargs -I % echo "`cat header.txt`\t%" > ${name}
    
    # Avoid intersecting with empty file
    if [[ `wc -l <${snps_file}` -ge 2 ]]; then
	    bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb >> ${name}
    fi
	"""
}

process split_into_samples {
    tag "${indiv_id}"
    conda params.conda

    input:
        tuple val(indiv_id), path(pval_file)
    
    output:
        path "*${suffix}"
    
    script:
    suffix = "sample_split.bed"
    """
    python3 $moduleDir/bin/split_into_samples.py ${pval_file} ${suffix}
    """
}


process annotate_variants {
    conda params.conda
    tag "${sample_id}"
    publishDir "${params.outdir}/by_sample"
    scratch true
    label "med_mem"

    input:
        tuple val(sample_id), path("pvals.bed"), path(peaks_file), path(footprint_file), path(hotspots_file)

    output:
        tuple val(sample_id), path(name)

    script:
    name = "${sample_id}.nonaggregated.bed"
    """
    process_file() {
        # \$1 - pval file, no header
        # \$2 - peaks file, no header
        # \$3 - output file

        bname=`basename \$2`

        if [[ "\${bname}" == *".empty" ]]; then
            # Nextflow needs an input file. Mock file was passed here.
            awk '{print "-"}' \$1 > \$3
        else
            if [[ "\${bname}" != *".bed.gz" ]]; then
                (bedmap --indicator \$1 \$2 > \$3) || true
            else
                (zgrep -v '#' \$2 | bedmap --indicator \$1 - > \$3) || true
            fi
        fi
    }

    tail -n+2 pvals.bed | sort-bed - > pval_f.sorted.bed
    head -1 pvals.bed > pvals.header.txt

    process_file \
        pval_f.sorted.bed \
        ${footprint_file} \
        footprints.txt

    process_file \
        pval_f.sorted.bed \
        ${peaks_file} \
        peaks.txt

    process_file \
        pval_f.sorted.bed \
        ${hotspots_file} \
        hotspots.txt

    cat pvals.header.txt pval_f.sorted.bed > sorted_pvals.bed

    # Doing it here to check if footprints and hotspots columns are present
    python3 $moduleDir/bin/add_annotations.py \
        --hotspots hotspots.txt \
        --peaks peaks.txt \
        --footprints footprints.txt \
        sorted_pvals.bed \
        ${name}
    """
}

workflow estimateBAD {
    take:
        snv_files
        prefix
    main:
        non_empty = snv_files
            | filter { it[1].countLines() >= params.min_snps_count }
        
        // At least one line except header
        out = apply_babachi(non_empty, prefix).intersect
            | filter { it[1].countLines() > 1 }
    emit:
        out
}

workflow combineFilesAndEstimateBAD {
    take:
        data
        prefix
    main:
        babachi_files = data
            | groupTuple()
            | collect_files

        out = estimateBAD(babachi_files, prefix) 
    emit:
        babachi_files
        out
}

workflow {
    // Estimate BAD and call 1-st round CAVs
    bads = Channel.of(params.states.tokenize(','))
    input_data = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(row.indiv_id, file(row.snps_file)))
        | filter { it[1].exists() }

    intersect_files = combineFilesAndEstimateBAD(input_data, 'iter1')[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFirstRound(intersect_files)

    // Reestimate BAD, and add excluded SNVs
    all_snps = estimateBAD(no_cavs_snps, 'final')
        | join(intersect_files, remainder: true)
        | map(it -> tuple(it[0], it[1] != null ? it[1] : file('empty'), it[2]))
        | add_cavs

    // Annotate with footprints and hotspots + aggregate by provided aggregation key
    out = calcPvalBinom(all_snps)
        | split_into_samples
        | flatten()
        | map(it -> tuple(it.name.replaceAll('.sample_split.bed', ''), it))
        | annotateWithFootprints
    
    aggregation(out)
}  

workflow annotateWithFootprints {
    take:
        pval_files
    main:
        annotations = Channel.fromPath(params.samples_file)
            | splitCsv(header: true, sep: '\t')
            | map(row -> tuple(
                    row.sample_id,
                    check_var(row['peaks_file_0.001fdr'], 'peaks'), 
                    check_var(row?.footprints_file, 'fp'),
                    check_var(row['hotspots_file_0.05fdr'], 'hs')
                    )
                )
        out = pval_files
            | join(annotations)
            | annotate_variants
    emit:
        out
}


workflow reannotateWithFootprints {
    Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))
        | annotateWithFootprints
}
