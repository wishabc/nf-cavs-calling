#!/usr/bin/env nextflow


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
    ${command} ${pval_files} \
        | cut -f1-6 \
        | sort-bed - \
        | uniq >> ${name}
    """
}


process ld_scores {
	conda params.conda
	tag "${prefix}"

	input:
		tuple val(chrom), path(snps_positions)
	
	output:
		path name
	
	script:
    prefix = "${chrom}:${snps_positions.simpleName}"
	name = "${prefix}.geno.ld"
    additional_params = chrom == 'all' ? "" : "--chr ${chrom}"
 	"""
    echo -e "chrom\tchromStart\tchromEnd" > variants.bed
    cat ${snps_positions} \
        | awk -v OFS='\t' '{ print \$1,\$2,\$3 }'  \
        | sort-bed - \
        | uniq >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
		--minDP 10 \
        --ld-window-bp ${params.ld_window_size_bp} \
        --bed variants.bed \
		--out ${prefix} \
        ${additional_params} \
	"""
}

process sort {
    conda params.conda
    label "high_mem"

    input:
        path ld_scores
    
    output:
        path name
    
    script:
    name = "${ld_scores.baseName}.sorted.bed"
    """
    tail -n+2 ${ld_scores} \
        | awk -v OFS='\t' '{print \$1,\$2-1,\$2,\$3,\$4,\$5}' \
        | sort-bed - > ${name}
    """
}

process intersect_with_tested_variants {
    conda params.conda
	tag "${prefix}"

    input:
        tuple path(ld_scores), path(variants_file)
    
    output:
        path name
    
    script:
    prefix = variants_file.simpleName
    name = "${prefix}.bed"
    """
    # Format: chr, start, end, end2, distance, start_es, end_es, sample, chr, start, end, end3 ...
    echo -e "#chr\tend1\tend2\tdistance\tes1\tes2\tsample_id\tn_indiv\tr_squared" > ${name}
    python3 $moduleDir/bin/find_neighbors.py ${variants_file} \
        | sort-bed - \
        | bedtools intersect -a stdin -b ${ld_scores} -wa -wb -sorted \
        | awk -v OFS='\t' '\$4==\$12 { print; }' \
        | cut -f-1,3-8,13- >> ${name}
    """
}


workflow annotateLD {
    take:
        data
        samples
    main:
        out = Channel.of(1..22)
            | map(it -> "chr${it}")
            | combine(data)
            | ld_scores
            | collectFile(
                keepHeader: true,
                sort: true,
                skip: 1,
                name: "ld_scores.geno.ld"
            )
            | sort
            | combine(samples)
            | intersect_with_tested_variants
            | collectFile(
                storeDir: params.outdir,
                keepHeader: true,
                skip: 1,
                name: "ld_scores.annotated_samples.geno.ld"
            )
    emit:
        out
}

workflow {
    by_sample = Channel.fromPath("${params.main_run_outdir}/by_sample/*.bed") 
    annotateLD(
        by_sample
            | collect(sort: true)
            | filter_tested_variants,
        by_sample
    )
}
