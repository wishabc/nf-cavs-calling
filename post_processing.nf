#!/usr/bin/env nextflow
include { filter_tested_variants } from "./main"
include { annotateLD } from "./ld_scores"
include { cavsMotifEnrichment } from "./motif_enrichment"


process random_sample {
    tag "seed:${step_start}-${step_start+params.samples_per_job}"
    conda params.conda
    label "sampling"

    input:
        tuple val(step_start), path(annotations_file), path(non_aggregated_file)

    output:
        path name

    script:
    name = "${step_start}.sampling.tsv"
    """
    python3 $moduleDir/bin/random_sample.py \
        ${non_aggregated_file} \
        ${annotations_file} \
        ${name} \
        --start ${step_start} \
        --step ${params.samples_per_job}
    """
}


// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"

process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda
    label "med_mem"

    input:
        tuple path(vcf), path(variants_file)

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

process extract_context {
    conda params.conda
    scratch true

    input:
        path variants
    output:
        path name

    script:
    name = "variants_context.bed"
    """
    echo -e "#chr\tstart\tend\tsequence" > ${name}
    cat ${variants} \
        | awk -v OFS='\t' '{ print \$1,\$2-${params.window},\$3+${params.window} }' \
        | uniq > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -bedOut \
        | awk -v OFS='\t' '{ print \$1,\$2+${params.window},\$3-${params.window},\$4 }' >> ${name}
    """
}

// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    publishDir params.outdir
    label "high_mem"

    input:
        path pval_file

    output:
        path name
        
    when:
        file(params.phenotypes_data, type: 'dir').exists()

    script:
    name = "phenotypes_ann.bed"
    """
    echo "Annotating"
    python3 $moduleDir/bin/annotate_with_phenotypes.py \
        ${params.phenotypes_data} \
        ${pval_file} \
        ${name}
    """
}

process merge_annotations {
    conda params.conda
    publishDir params.outdir
    scratch true

    input:
        path unique_snps
        path context
        path mutation_rates
    
    output:
        path name
    
    script:
    name = "cavs.annotations.bed.gz"
    """
    python3 $moduleDir/bin/merge_annotations.py \
        ${unique_snps} \
        ${context} \
        ${mutation_rates} \
        tmp.bed
    
    head -1 tmp.bed > res.bed
    sort-bed tmp.bed >> res.bed
    bgzip -c res.bed > ${name}
    """
}

workflow sampleVariants {
    take:
        data
    main:
        params.sampling_count = 100
        params.samples_per_job = 100
        total_count = params.sampling_count * params.samples_per_job 
        out = Channel.of(1..params.sampling_count) 
            | map(it -> it * params.samples_per_job)
            | combine(data)
            | random_sample
            | collectFile(
                storeDir: params.outdir,
                name: "subsampled.n${total_count}.tsv",
                keepHeader: true,
                skip: 1
            )
    emit:
        out
        
}


workflow mutationRates {
    take:
        data
    main:
        out = Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
            | combine(data)
            | process_mutation_rates
            | collectFile(
                name: "mut_rates.annotation.bed",
                skip: 1,
                keepHeader: true
            )
    emit:
        out
}

workflow {
    // sample_wise_pvals = Channel.fromPath("${params.raw_pvals_dir}/*.bed")

    // data = sample_wise_pvals
    //     | collect(sort: true)
    //     | filter_tested_variants
    data = Channel.fromPath("${params.data}")
    //annotateLD(sample_wise_pvals, data)

    data | (annotate_with_phenotypes & cavsMotifEnrichment)

    merge_annotations(data, extract_context(data), mutationRates(data))

}

workflow sample {
    nonagr_files = Channel.fromPath(params.nonagr_pvals)
    Channel.fromPath("${params.outdir}/cavs.annotations.bed.gz")
        | combine(nonagr_files)
        | sampleVariants
}