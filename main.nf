#!/usr/bin/env nextflow
include { estimateBadByIndiv; estimateBad; intersectWithBadmap } from "./bad_estimation"
//include { extractAndFilter } from "./extract_and_filter"
include { callCavsFromVcfsBinom; calcPvalBinom; calcPvalNegbin; fitNegBinom } from "./pval_calc"

process concat_files {
  input:
  tupple val(id) path(x)
  val out_file
  output:
  file out_file
  script:
  """
  < $x zcat > ${out_file}
  """
}

workflow {
    intersect_map = estimateBadByIndiv()
    no_cavs_snps = callCavsFromVcfsBinom(intersect_map)
    new_badmap = estimateBad(no_cavs_snps, 'nocavs_badmap')
    new_badmap_join = intersect_map.join(new_badmap)
    new_intersect_map = intersectWithBadmap(new_badmap_join, 'nocavs_intersect')
    switch (params.strategy) {
        case 'binom':
            calcPvalBinom(new_intersect_map, 'nocavs_')
            break
        case 'negbin':
            badmaps = concat_files(new_intersect_map, 'badmaps.tsv').collectFile()
            badmaps.view()
            //weights_files = fitNegBinom(badmaps)
            //calcPvalNegbin(new_intersect_map, weights_files, 'nocavs_')
            break
        default:
            println("Wrong strategy provided. ${strategy} not in ('binom', 'negbin')")
    }
    
}