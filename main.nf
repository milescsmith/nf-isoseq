#!/usr/bin/env nextflow
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


Channel
    .from( "$params.raw_dir/*.ccs.subreads.bam" , checkIfExists:true, type: 'glob' )
    .set{raw_subreads}

Channel
    .fromPath(params.barcode)
    .splitCsv(header:false)
    .map{ row -> row[0].split("-")[0] }
    .unique()
    .set{sample_channel}


process ccs_calling {
    tag "CCS calling"
    publishDir '${params.css_dir}', mode: 'copy', pattern: '*.ccs.bam', overwrite: true
    publishDir '${params.css_dir}', mode: 'copy', pattern: "*.log", overwrite: true

    input:
        file subreads from raw_subreads

    output:
        file "*.ccs.bam" into ccs_channel
        file "css.log" into css_log

    script:
        """
        ccs \
            --min-rq 0.9 \
            --log-level DEBUG \
            --log-file ccs.log \
            --num-threads ${task.cpus} \
            ${subreads} \
            ${params.runid}.ccs.bam
        """

}

process demultiplex {
    tag "Demultiplexing"
    publishDir '${params.demux_dir}', mode: 'copy', pattern: '*.bam', overwrite: true
    publishDir '${params.demux_dir}', mode: 'copy', pattern: '*.log', overwrite: true

    input:
        called_ccs from ccs_channel
    
    output:
        file "*.bam" into demuxed_channel
        file "*.log" into demux_log

    script:
        """
        lima \
            --isoseq \
            --dump-clips \
            --peek-guess \
            --num-threads ${task.cpus} \
            --log-level INFO \
            --log-file lima.demux.log \
            ${called_ccs} \
            barcode.primers.fasta \
            ${demux_dir}/demuxed.bam
        """
}

process refine {
    tag "Refining $sample"
    publishDir '${params.refine_dir}', mode: 'copy', pattern: '*.flnc.bam', overwrite: true

    input:
        val sample from sample_channel
        demuxed_files from demuxed_channel
    
    output:

    script:
    """
    isoseq3 \
        refine \
        --num-threads ${tasks.cpus} \
        --log-file ${sample}.log \
        --log-level INFO \
        --verbose \
        ${demuxed_dir}/demuxed.${sample}.bam \
        primers.fasta \
        ${refined_dir}/refined.${sample}.flnc.bam
    """
}

process cluster-polish {
    tag "Polishing $sample"
    publishDir '${params.polished_dir}', mode: 'copy', pattern: '*.bam', overwrite: true

    input:
        val sample from sample_channel
    
    output:

    script:
    """
    isoseq3 \
        cluster \
        --use-qvs \
        --num-threads ${task.cpu} \
        --log-file ${sample}.log \
        --log-level INFO \
        --verbose \
        ${refined_dir}/refined.${sample}.flnc.bam \
        ${polished_dir}/polished.${sample}.bam
    """
}

process mapping {
    tag "Mapping $sample"
    publishDir '${params.mapped_dir}', mode: 'copy', pattern: '*.sam', overwrite: true
    publishDir '${params.mapped_dir}', mode: 'copy', pattern: '*.log', overwrite: true

    input:
        val sample from sample_channel

    output:

    script:
    """
    gmap \
        -D ${params.reference}
        -d ${species} \
        -f samse \
        -n 0 \
        -t ${task.cpus} \
        --cross-species \
        --max-intronlength-ends 200000 \
        -z sense_force \
        ${polished_dir}/polished.${sample}.hq.fasta \
        > ${mapped_dir}/${sample}_hq_isoforms.mapped.sam \
        2> ${mapped_dir}/${sample}_hq_isoforms.log
    """
}

process sort {
    tag "Sorting"
    publishDir '${params.sorted_dir}', mode: 'copy', pattern: '*.bam', overwrite: true

    input:
        val sample from sample_channel

    output:
    
    script:
    """
    samtools \
        sort \
        -O BAM \
        ${mapped_dir}/${sample}_hq_isoforms.mapped.sam \
        -o ${sorted_dir}/${sample}_hq_isoforms.mapped.sorted.bam 
    """
}

process transcompress {
    tag "Transcompressing $sample"
    publishDir '${params.transcompressed_dir}', mode: 'copy', pattern: '*.fastq.gz', overwrite: true

    input:
        val sample from sample_channel

    output:

    script:
    """
    bgzip \
        --decompress ${polished_dir}/${sample}.fasta.gz && \
    bgzip \
        --index ${polished_dir}/${sample}.fasta
    """
}

process filter {
    tag "Filtering $sample"
    publishDir '${params.filtered_dir}', mode: 'copy', pattern: '*.bam', overwrite: true

    input:
        val sample from sample_channel

    output:

    script:
    """
    filter_sam \
        --fastq ${polished_dir}/${sample}.fastq.gz \
        --sam ${sorted_dir}/${sample}_hq_isoforms.mapped.sorted.bam \
        --prefix ${filtered_dir}/${sample}_
    """
}

process collapse {
    tag "Collapsing"
    publishDir '${params.collapsed_dir}', mode: 'copy', pattern: '*.ccs.bam', overwrite: true

    input:
        val sample from sample_channel

    output:
    
    script:
    """
    collapse_isoforms_by_sam.py \
        --input ${polished_dir}${sample}.fastq \
        --sam ${sorted_dir}/${sample}_filtered.sam \
        --dun-merge-5-shorter \
        --prefix ${collapsed_dir}/${sample}
    """
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}