#!/usr/bin/env nextflow

params.mode = filter
params.samples_file = '/home/sabramov/nf-babachi/metadata+ag_number.tsv'
params.outdir = '/home/sabramov/nf-babachi/babachi_all_states/'
params.vcf_file = "/home/sabramov/nf-babachi/allele_counts.fixed.vcf.gz"
params.states = "1,1.5,2,2.5,3,4,5,6"
params.prior = "uniform"
params.filtered_vcfs = ''

// DO NOT EDIT BELOW

process extract_indiv_vcfs {
    tag "Extracting ${indiv_id}"
    publishDir params.outdir + 'indiv_vcfs'

    cpus 5

    input:
	    tuple val(indiv_id), val(agg_numbers)
    output:
        tuple val(indiv_id), path("${indiv_id}.snps.bed")
    script:
    """
    bcftools view -s ${agg_numbers} --threads ${task.cpus} ${params.vcf_file} > ${indiv_id}.vcf
    """
}

process filter_indiv_vcfs {
    tag "Filtering ${indiv_id}"
    publishDir params.outdir + 'filtered_indiv_vcfs'

    input:
	    tuple val(indiv_id), path(indiv_vcf)
    output:
        tuple val(indiv_id), path(get_filtered_file_by_indiv_id(indiv_id))
    script:
    name = get_filtered_file_by_indiv_id(indiv_id)
    """
    babachi filter ${indiv_vcf} -O ${name}
    rm ${indiv_vcf}
    """
}

def get_filtered_file_by_indiv_id(indiv_id) {
    "${indiv_id}.snps.bed"
}

process apply_babachi {
	cpus 2
    publishDir params.outdir + '/badmaps'

	input:
		tuple val(indiv_id), path(snp_file)
	output:
		tuple val(indiv_id), path(snp_file), path("${indiv_id}.bad.bed")

	script:
	"""
    babachi ${snps_file} -O ${indiv_id}.bad.bed --visualize -z -e png \
	 -j ${task.cpus} -p ${prior} -s ${states}
	"""
}


process intersect_with_snps {

    publishDir params.outdir + '/snp_annotation'
	input:
		tuple val(indiv_id), path(snps_file), path(badmap_file)
    
    output:
        tuple val(indiv_id), path("${indiv_id}.intersect.bed")

	script:
	"""
	bedtools intersect -a ${snps_file} -b ${badmap_file} -wa -wb | awk \
    '{for (i=1;i<=7;i+=1) print \$i, print \$11}' > ${indiv_id}.intersect.bed
	"""
}


//indiv_badmap_intersect
//    .map{ it -> it[1] }
//    .collectFile(name: 'bad_annotations_files.txt', newLine: true)
//    .set{bad_annotations}


// process collect_stats_for_neg_bin {
//     input:
//         path bad_annotations
//     output:
//         path "nb_fit" into nb_fit
    
//     script:
//     """
//     python3 collect_nb_stats.py ${bad_annotations} ${nb_stats_dir}
//     """

// }

// process annotate_with_db {
//     input:
//         tuple val(indiv_id), path(snps_file) from indiv_snps_file

//     output:
//         tuple val(indiv_id), path('')

//     script:
//     """
//     echo ANNOTATE WITH COSMIC
//     """
// }

// process add_read_counts {
//     input:
//         tuple val(indiv_id), path(badmap_intersect_file) from indiv_badmap_intersect

//     script:
//     """
//     echo ADD READ COUNTS ${indiv_id}
//     """
// }

// process calculate_bin_pvalue {
//     input:
//         tuple val(indiv_id), path(badmap_intersect_file) from indiv_badmap_intersect
//     output:
//         tuple val(indiv_id), path("${indiv_id}.pvalue") into indiv_pvalues

//     script:
//     """
//     echo Calc bin pval ${indiv_id}
//     """
// }

// process calculate_nbin_pvalue {
//     input:
//         tuple val(indiv_id), path(pvalue_file) from indiv_pvalues
//         path nb_fit
//     output:
//         tuple val(indiv_id), path("${indiv_id}.negbin.pvalue") into indiv_pvalues

//     script:
//     """
//     echo Calc negbin pval ${indiv_id} ${nb_fit}
//     """
// }

workflow extract_and_filter {
    main:
        sample_ag_merge = Channel
                .fromPath(params.samples_file)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id, row.ag_number) }
                .groupTuple(by:0)
                .map{ it -> tuple(it[0], it[1].join(",")) }
                .last()
        println(sample_ag_merge)
        extract_indiv_vcfs(sample_ag_merge) | filter_indiv_vcfs 
    emit:
        filter_indiv_vcfs.out
}

workflow {
    if ( params.filtered_vcfs != '' ) 
        extracted_vcfs = Channel.fromPath(params.samples_file)
                .splitCsv(header:true, sep:'\t')
                .map{ row -> tuple(row.indiv_id, path(get_filtered_file_by_indiv_id(row.indiv_id))) }
    else   
        extracted_vcfs = extract_and_filter()
    
    //apply_babachi(extracted_vcfs) | intersect_with_snps
}