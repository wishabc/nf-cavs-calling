#!/usr/bin/env nextflow

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
/*
workflow calcPval {
    main:
        extracted_vcfs = Channel.fromPath(params.samples_file)
            .splitCsv(header:true, sep:'\t')
            .map{ row -> tuple(row.indiv_id,
                path(params.filtered_vcfs + '/' + get_filtered_file_by_indiv_id(row.indiv_id))) }

        apply_babachi(extracted_vcfs) | intersect_with_snps
    emit:
        intersect_with_snps.out
}

workflow {
    calcPval()
}
*/