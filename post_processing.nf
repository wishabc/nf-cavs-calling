// DEFUNC

// Distance to the nearest peak for tested variants
// Common variants

process distance_to_dhs {

    input:
        tuple path(variants), val(dhs)
    
    
    script:
    """
    cut -f1-3 ${params.dhs_masterlist} \
        | closest-features \
            --closest \
            --dist \
            --delim '\t' \
            ${variants} \
            - > ${name}
    """
}