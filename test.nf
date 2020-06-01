echo true

// params.barcode = "/mnt/scratch/guth-aci/pacbio/OMRF_pool/barcode_oligos.csv"

extract_bc = { item -> 
    item =~ /bc(\d+)\-\[FR]/
}

params.input = "stuff"

Channel
    .from( params.input )
    .set{ raw_subreads_1 }

Channel
    .of( 1..params.chunks )
    .set{ chunk_list }

process printForward {

    input: 
        val x from raw_subreads_1
        each y from 1..params.chunks

    script:
        """
        echo '${x} ${y}'
        """

}