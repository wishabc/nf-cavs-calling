#!/usr/bin/env nextflow
include { filter_tested_variants } from "./main"
include { annotateLD } from "./ld_scores"
include { cavsMotifEnrichment } from "./motif_enrichment"
include { differentialCavs } from "./differential_cavs"


// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"
params.phenotypes_data = "/home/sabramov/phenotypes_data"


process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda
    label "med_mem"

    input:
        tuple  path(vcf), path(variants_file)

    output:
        path name

    script:
    name = "${vcf.simpleName}.bed"
    """
    echo -e "#chr\tstart\tend\tID\tref\talt\tchr\tstart_mr\tend_mr\tref_mr\talt_mr\tmut_rates_roulette\tmut_rates_gnomad" > tmp.bed
    
    bcftools query -f"%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/MR\t%INFO/MG\n" \
        ${vcf} | awk '{print "chr"\$0}' | bedtools intersect \
        -a ${variants_file} -b stdin -sorted -wa -wb >> tmp.bed
    
    python3 $moduleDir/bin/filter_variants.py tmp.bed ${name}
    """
}

process merge_and_sort {
    publishDir "${params.outdir}"
    conda params.conda
    scratch true

    input:
        path bed_files
    
    output:
        path name

    script:
    name = "mut_rates.annotation.bed"
    """
    for file in ${bed_files}; do
        awk 'NR>1' \$file >> tmp.bed
    done
    head -1 ${bed_files[0]} > ${name}
    sort-bed tmp.bed >> ${name}
    """
}

process extract_context {
    conda params.conda
    scratch true
    publishDir "${params.outdir}"

    input:
        path variants
    output:
        path name

    script:
    name = "variants_context.bed"
    """
    cat ${variants} \
        | awk -v OFS='\t' '{ print \$1,\$2-${params.window},\$3+${params.window} }' \
        | uniq > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -bedOut \
        | awk -v OFS='\t' '{ print \$1,\$2+${params.window},\$3-${params.window},\$4 }' > ${name}
    """
}

// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path pval_file

    output:
        path name

    script:
    name = "phenotypes_ann.bed"
    """
    python3 $moduleDir/bin/annotate_with_phenotypes.py ${params.phenotypes_data} ${pval_file} ${name}
    """
}


workflow mutationRates {
    take:
        data
    main:
        out = Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
            | combine(data)
            | process_mutation_rates
            | collect(sort: true)
            | merge_and_sort
    emit:
        out
}

workflow {
    Channel.fromPath(params.nonagr_pvals)
        differentialCavs

    sample_wise_pvals = Channel.fromPath("${params.raw_pvals_dir}/*.bed")

    data = sample_wise_pvals
        | collect(sort: true)
        | filter_tested_variants
    
    annotateLD(sample_wise_pvals, data)
    
    cavsMotifEnrichment()

    // merge_results(
        extract_context(data)
        mutationRates(data)
        annotate_with_phenotypes(data)
    // )

}
