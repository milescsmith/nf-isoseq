/*
 * Copyright (c) 2020, Oklahoma Medical Research Foundation (OMRF).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 *
 */

// Nextflow pipeline for processing PacBio IsoSeq runs
// Author: Miles Smith <miles-smith@omrf.org>
// Date: 2020/05/26
// Version: 0.1.0

// File locations and program parameters

log.info """\
Isoseq3-NF Pipeline
===================
params.base             = "/s/guth-aci/isoseq"
params.reference        = "/Volumes/guth_aci_informatics/references/genomic/homo_sapiens/indices/gmap/gencode_v32/homo_sapiens/"
params.barcodes         = "${params.base}/barcode_oligos.csv"
"""

Channel
    .fromFile( params.ccs , checkExists:true)
    .into{raw_subreads}

process ccs_calling{
    tag "CCS calling"
    publishDir "${params.base}/02_ccs", mode: "copy", pattern: "*.ccs.bam", overwrite: true

    input:
    set file(subreads) from raw_subreads

    output:
    file "*.ccs.bam" into ccs_channel

    script:
    """
    ccs --min-rq 0.9 --log-level DEBUG --log-file ccs.log --num-threads ${task.cpus} ${subreads} ${runid}.ccs.bam
    """


}

process demultiplex{

    script:
    """
    lima --isoseq --dump-clips --peek-guess --num-threads 36 --log-level INFO --log-file lima.demux.log $PROJECT_FOLDER/02_CSS/{runid}.ccs.bam barcode.primers.fasta $PROJECT_FOLDER/03_demultiplexed/demuxed.bam
    """
}

process refine{

    script:
    """
    isoseq3 refine --num-threads 16 --log-file ${bc_array[$i]}.log --log-level INFO --verbose $PROJECT_FOLDER/03_demultiplexed/demuxed.${bc_array[$i]}.bam primers.fasta $PROJECT_FOLDER/04_refined/refined.${bc_array[$i]}.flnc.bam
    """
}

process cluster-polish{

    script:
    """
    isoseq3 cluster --use-qvs --num-threads 16 --log-file ${bc_array[$i]}.log --log-level INFO --verbose $PROJECT_FOLDER/04_refined/refined.${bc_array[$i]}.flnc.bam $PROJECT_FOLDER/05_polished/polished.${bc_array[$i]}.bam
    """
}

process mapping{

    script:
    """
    gmap -D /Volumes/guth_aci_informatics/references/genomic/homo_sapiens/indices/gmap/gencode_v32/homo_sapiens/ -d homo_sapiens -f samse -n 0 -t 16 --cross-species --max-intronlength-ends 200000 -z sense_force $PROJECT_FOLDER/05_polished/polished.${bc_array[$i]}.hq.fasta > $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.fasta.sam 2> $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.log
    """
}

process sort{

    script:
    """
    samtools sort -O BAM $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.mapped.bam -o $PROJECT_FOLDER/07_sorted/${bc_array[$i]}_hq_isoforms.mapped.sorted.bam
    """
}

process transcompress{

    script:
    """
    bgzip --decompress $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq.gz && bgzip --index $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq &&
    """
}
process filter{

    script:
    """
    filter_sam --fastq $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq.gz --sam $PROJECT_FOLDER/06_sorted/${bc_array[$i]}_hq_isoforms.mapped.sorted.bam --prefix $PROJECT_FOLDER/07_filtered/${bc_array[$i]}_
    """
}

process collapse{

    script:
    """
    collapse_isoforms_by_sam.py --input $PROJECT_FOLDER/06_polished/${bc_array[$i]}.fastq --sam $PROJECT_FOLDER/07_sorted/${bc_array[$i]}_filtered.sam --dun-merge-5-shorter --prefix $PROJECT_FOLDER/08_collapsed/${bc_array[$i]}
    """
}