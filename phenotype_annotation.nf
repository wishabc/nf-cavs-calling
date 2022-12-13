params.phenotypes_data = "/home/sabramov/phenotypes_data"


process annotate_with_phenotypes {
    conda params.conda
    tag "${sample_id}"
    publishDir "${params.outdir}/phenotypes"

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

// TODO wrap in apptainer
params.ldsc_conda = "/home/sabramov/miniconda3/envs/ldsc"

process find_ld {
    publishDir "${params.outdir}/l2"
    tag "${phen_name}"
    conda params.conda

    input:
        tuple val(phen_id), val(phen_name), path(sumstats_file), val(ld_prefix), path("ld_files/*")
    
    output:
        tuple val(phen_id), val(phen_name), path(sumstats_file), val(ld_prefix), path("ld_files/*")
    
    script:
    """
    /home/sabramov/projects/ENCODE4/ldsc/ldsc.py \
        --print-snps /net/seq/data2/projects/sabramov/LDSC/UKBB_hm3.snps.tsv \
        --ld-wind-cm 1.0 \
        --out ld_files/${ld_prefix} \
        --bfile /home/sabramov/LDSC/plink_files/1000G.EUR.hg38 \
        --annot ld_files/${ld_prefix} \
        --l2
    """
}

process run_ldsc {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc"
    tag "${phen_name}"

    input:
        tuple val(phen_id), val(phen_name) path(sumstats_file), val(ld_prefix), path("ld_files/*"), val(frq_prefix), path("frqfiles/*")
    
    output:
        tuple val(phen_id), val(phen_name), path("phen_results/${name}*")

    script:
    name = "${phen_id}.result"
    """

    mkdir phen_results
    /home/sabramov/projects/ENCODE4/ldsc/ldsc.py \
        --h2 ${phenotype_sumstats} \
        --ref-ld-chr ld_files/${ld_prefix} \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out phen_results/${name}
    """
}

workflow LDSC {
    params.ann_path = '/net/seq/data2/projects/sabramov/LDSC/test_intersection/baselineLD.'
    annotations = Channel.of(file(params.ann_path))
        .map(it -> tuple(it.name, file("${it}*")))

    phens = Channel.fromPath("/net/seq/data2/projects/sabramov/LDSC/UKBB.phenotypes.test.tsv")
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.phen_id, row.phen_name, file(row.sumstats_file)))
    
    params.frqfiles = "/home/sabramov/LDSC/plink_files/1000G"
    params.weights = "/home/sabramov/LDSC/weights/weights."
    ld_data = find_ld(phens.combine(annotations))
    frqs = Channel.of(file(params.frqfiles))
        .map(it -> tuple(it.name, file("${it}*.frq")))
    run_ldsc(ld_data.combine(frqs))
}


workflow annotateWithPheno {
    params.pval_file_dir = "/net/seq/data2/projects/sabramov/ENCODE4/cav-calling/babachi_1.5_common_final/all_aggregations/output/final.ag_files_binom.all"
    pvals = Channel.fromPath("${params.pval_file_dir}/*.bed")
        .map(it -> file(it))
    annotate_with_phenotypes(pvals)
    
}