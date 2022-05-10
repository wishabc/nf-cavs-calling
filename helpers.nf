def get_file_by_indiv_id(ind, file_type='vcf') {
    switch (file_type) {
        case "vcf":
            return "${ind}.vcf.gz"
        case "filter":
            return "${ind}.snps.bed"
        case "badmap":
            return "${ind}.bad.bed"
        case "intersect":
            return "${intersect}.intersect.bed"
        default:
            return "default"
    }

}

def get_id_by_sample(ind, agg_id) {
    return "${ind}@${agg_id}"
}