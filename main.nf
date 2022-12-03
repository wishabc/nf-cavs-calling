#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad } from "./bad_estimation"
include { callCavsFromVcfsBinom; calcPvalBinom; addImputedCavs; aggregate_pvals } from "./pval_calc"
include { motifEnrichment } from "./motif_enrichment"


params.conda = "$moduleDir/environment.yml"
params.sample_pvals_dir = "$launchDir/${params.outdir}/sample_pvals"
params.footprints_master = ""


process sort_and_gzip {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path inp

    output:
        tuple path(name), path("${name}.tbi")

    script:
    name = "${inp.simpleName}.sorted.bed.gz"
    """
    cat ${inp} | grep -v '^#' | sort-bed - | bgzip -c > ${name}
    tabix ${name}
    """
}


process split_into_samples {
    tag "${indiv_id}"
    conda params.conda
    publishDir "${params.outdir}/sample_pvals"

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
    publishDir "${params.outdir}/annotations"
    scratch true

    input:
        tuple val(sample_id), path(pval_file), path(hotspots_file), val(footprint_file)

    output:
        tuple val(sample_id), path(name)

    script:
    name = "${sample_id}.fp_annotation.bed"
    footprint_f = footprint_file ? footprint_file : "empty.bed"
    """
    if [[ "${footprint_f}" == "empty.bed" ]]; then
        touch ${footprint_f}
    fi
    sort-bed ${pval_file} > pval_f.bed

    bedmap --header \
        --indicator pval_f.bed \
        ${footprint_f} >> footprints.txt
    
    bedmap --header \
        --indicator pval_f.bed \
        ${hotspots_file} >> hotspots.txt

    echo -e "`head -1 ${pval_file}`\tfootprints\thotspots\n" > ${name}
    paste pval_f.bed footprints.txt hotspots.txt >> ${name}
    """
}


workflow aggregation {
    take:
        sample_split_pvals
    main:
        agg_key = params.aggregation_key ? params.aggregation_key : "all"
        store_dir = "${params.outdir}/pvals_nonaggregated.${agg_key}"
        if (agg_key != 'all') {
            sample_cl_correspondence = Channel.fromPath(params.samples_file)
                    .splitCsv(header:true, sep:'\t')
                    .map(row -> tuple(row.ag_id, row[params.aggregation_key]))
            pvals = sample_split_pvals
                .join(sample_cl_correspondence)
                .filter(it -> !it[2].isEmpty())
                .collectFile(keepHeader: true,
                 skip: 1, storeDir: store_dir) { item -> [ "${item[2]}.bed", item[1].text + '\n' ]}
                .map(it -> tuple(it.simpleName, it))
        } else {
            pvals = sample_split_pvals.map(it -> it[1])
                .collectFile(name: 'all_pvals.bed',
                storeDir: store_dir,
                keepHeader: true, skip: 1)
            .map(it -> tuple('all', it))
        }
        footprints = Channel.fromPath(params.footprints_master)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.ag_id, file(row.footprint_path)))
        ann_pvals = pvals.take(1)
        //ann_pvals = annotateWithFootprints(pvals, footprints)
        out = aggregate_pvals(ann_pvals, "binom.${agg_key}", 'final.') // | motifEnrichment
    emit:
        out[0]
}

workflow annotateWithFootprints {
    take:
        pval_files
        footprints
    main:
        data = pval_files.join(footprints, remainder: true)
        annotations = annotate_variants(data)
    emit:
        annotations
}


workflow withExistingFootprints {
    sample_pvals = Channel.fromPath("${params.sample_pvals_dir}/*.bed")
        .map(it -> tuple(file(it).simpleName, file(it)))
    hotspots = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.ag_id, file(row.hotspots_file)))
    d = sample_pvals.join(hotspots)

    footprints = Channel.fromPath(params.footprints_master)
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.ag_id, file(row.footprint_path)))
    annotateWithFootprints(d, footprints)
}

workflow aggregatePvals {
    sample_pvals = Channel.fromPath("${params.sample_pvals_dir}/*.bed")
        .map(it -> tuple(file(it).simpleName, file(it)))
    aggregation(sample_pvals)
}


workflow {
    // Estimate BAD and call 1-st round CAVs
    iter1_prefix = 'iter1.'

    bads = Channel.of(params.states.split(','))
    filtered_vcfs_and_intersect = estimateBadByIndiv(iter1_prefix)
    filtered_vcfs = filtered_vcfs_and_intersect[0]
    intersect_files = filtered_vcfs_and_intersect[1]
    // Calculate P-value + exclude 1-st round CAVs 
    no_cavs_snps = callCavsFromVcfsBinom(intersect_files, iter1_prefix)


    iter2_prefix = 'final.'
    // Reestimate BAD, and add excluded SNVs
    iter2_intersections = estimateBad(no_cavs_snps, iter2_prefix)
    imputed_cavs = addImputedCavs(iter2_intersections.join(intersect_files))
    binom_p = calcPvalBinom(imputed_cavs, iter2_prefix)[0]
    sample_split_pvals = split_into_samples(binom_p).flatten()
        .map(it -> tuple(it.simpleName, it))
    aggregation(sample_split_pvals)
    }