#!/usr/bin/env nextflow
include { callCavsFirstRound; calcPvalBinom; add_cavs } from "./pval_calc"
include { aggregation } from "./aggregation"

params.conda = "$moduleDir/environment.yml"


def set_key_for_group_tuple(ch) {
  ch.groupTuple()
    | map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
    | transpose()
}

def check_var(var, prefix) {
    if (var) {
        if (file(var).exists()) {
            return file(var)
        }
    } 
    return file("${prefix}.empty") 
}

process filter_tested_variants {
    conda params.conda
    scratch true

    input:
        path pval_files

    output:
        path name

    script:
    // Expected all files to be in the same format
    command = pval_files[0].extension == 'gz' ? 'zcat' : 'cat'
    name = pval_files.size() > 1 ? "unique_variants.bed" : "${pval_files[0].simpleName}.bed"
    """
    head -1 ${pval_files[0]} | cut -f1-6 > ${name}
    
    ${command} ${pval_files} \
        | awk -v OFS='\t' -v col='is_tested' \
            'NR==1 {
                for(i=1;i<=NF;i++){
                    if (\$i==col){
                        c=i;
                        break
                    }
                }
            }
            (\$c == "True") { print }' \
        | cut -f1-6 \
        | sort-bed - \
        | uniq >> ${name}
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

    babachi snps.common.bed \
        -O ${badmap_file} \
        -j ${task.cpus} \
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

    input:
        tuple val(sample_id), path(pval_file), path(hotspots_file), path(footprint_file)

    output:
        tuple val(sample_id), path(name)

    script:
    name = "${sample_id}.nonaggregated.bed"
    """
    sort-bed ${pval_file} > pval_f.bed

    # Add footprints
    if [[ "${footprint_file.name}" == "fp.empty" ]]; then
        cat pval_f.bed | awk '{ print "-" }' > footprints.txt
    else
        bedmap --header \
            --indicator pval_f.bed \
            ${footprint_file} > footprints.txt
    fi

    if [[ "${hotspots_file.name}" == "hp.empty" ]]; then
        cat pval_f.bed | awk '{ print "-" }' > hotspots.txt
    else
        bedmap --header \
            --indicator pval_f.bed \
            ${hotspots_file} > hotspots.txt
    fi

    echo -e "`head -1 ${pval_file}`\tfootprints\thotspots" > ${name}
    paste pval_f.bed footprints.txt hotspots.txt >> ${name}
    """
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
        | awk '\$7+\$8 >= ${params.initial_filter} {print;}' \
        | sort-bed - >> ${name}
    """
}

process estimate_mse {
    conda params.conda
    scratch true
    publishDir params.outdir
    label "high_mem"

    input:
        path all_pvals
    
    output:
        path name
    
    script:
    name = "mse_estimates.tsv"
    """
    python3 $moduleDir/bin/estimate_mse.py ${all_pvals} ${name}
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
    iter1_prefix = 'iter1'

    bads = Channel.of(params.states.tokenize(','))
    input_data = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(row.indiv_id, file(row.snps_file)))

    intersect_files = combineFilesAndEstimateBAD(input_data, iter1_prefix)[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFirstRound(intersect_files, iter1_prefix)

    iter2_prefix = 'final'
    // Reestimate BAD, and add excluded SNVs
    all_snps = estimateBAD(no_cavs_snps, iter2_prefix)
        | join(intersect_files, remainder: true)
        | map(it -> tuple(it[0], it[1] != null ? it[1] : file('empty'), it[2]))
        | add_cavs

    // Annotate with footprints and hotspots + aggregate by provided aggregation key
    out = calcPvalBinom(all_snps, iter2_prefix)
        | split_into_samples
        | flatten()
        | map(it -> tuple(it.name.replaceAll('.sample_split.bed', ''), it))
        | annotateWithFootprints

    mse = out
        | map(it -> it[1])
        | collectFile(name: "all_pvals.bed", skip: 1, keepHeader: true)
        | estimate_mse
    
    aggregation(out, mse)
}  

workflow annotateWithFootprints {
    take:
        pval_files
    main:
        annotations = Channel.fromPath(params.samples_file)
            | splitCsv(header: true, sep: '\t')
            | map(row -> tuple(row.ag_id,
                    check_var(row?.hotspot_peaks_point1per, 'hp'), 
                    check_var(row?.footprint, 'fp')
                    )
                )
        out = pval_files
            | join(annotations)
            | annotate_variants
    emit:
        out
}


workflow tmp {
    Channel.fromPath("${params.raw_pvals_dir}/*.bed")
        | map(it -> tuple(it.name.replaceAll('.nonaggregated.bed', ""), it))
        | annotateWithFootprints
}
