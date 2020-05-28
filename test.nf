echo true

params.barcode = "/mnt/scratch/guth-aci/pacbio/OMRF_pool/barcode_oligos.csv"

extract_bc = { item -> 
    item =~ /bc(\d+)\-\[FR]/
}

Channel
    .fromPath(params.barcode)
    .splitCsv(header:false)
    .map{ row -> row[0].split("-")[0] }
    .unique()
    .set{csv_channel}

process printForward {

    input: 
        val x from csv_channel

    script:
        """
        echo '${x}'
        """

}