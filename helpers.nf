def get_stats_dir() {
    stats_dir = "${params.outdir}/stats"
}
def get_file_by_indiv_id(ind, file_type) {
    switch (file_type) {
        case "vcf":
            return "${ind}.vcf.gz"
        case "badmap":
            return "${ind}.bad.bed"
        case "intersect":
            return "${ind}.intersect.bed"
        case "nocavs_intersect":
            return "${ind}.nocavs.intersect.bed"
        case "pvalue":
            return "${ind}.pvalue.bed"
        case "aggregation":
            return "${ind}.aggregation.bed"
        case "nocavs":
            return "${ind}.nocavs.bed"
        case "add_cavs":
            return "${ind}.added_cavs.intersect.bed"
        case "*":
            return "${ind}*"
        default:
            return "default"
    }

}

def get_id_by_sample(ind, agg_id) {
    return "${ind}@${agg_id}"
}