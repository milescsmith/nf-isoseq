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

params.input = "${params.raw}/*.subreads.bam"

Channel
    .fromPath( params.input , checkIfExists: true )
    .into{ raw_subreads_1; raw_subreads_2 }

Channel
    .value()
    .fromPath( params.barcodes )
    .splitCsv( header:false )
    .map{ row -> row[0].split("-")[0] }
    .unique()
    .set{ sample_name_ch }

Channel
    .fromPath( params.barcodes )
    .into{ barcodes_ch; refine_barcodes_ch }

extract_bc = { item -> 
    item =~ /bc(\d+)\-\[FR]/
}

// change most of these "bioconda::[package]" to pb.yml
process ccs_indexing {
    conda "bioconda::pbbam"

    tag "CCS indexing"
    // publishDir "${params.raw}", mode: "copy", pattern : "*.bam.pbi", overwrite: false
    
    input:
        file(subreads) from raw_subreads_1

    output:
        file "*.bam.pbi" into raw_index_ch

    script:
        """
        pbindex ${subreads}
        """
}

process ccs_calling {
    // Need to rewrite so as to chunk and parallelize https://github.com/PacificBiosciences/ccs#how-can-I-parallelize-on-multiple-servers

    // Okay, now make chunking optional
    conda "bioconda::pbccs"
    
    cpus 8
    maxForks 8

    tag "CCS calling"
    // publishDir "${params.ccs}", mode: "copy", pattern: "*.ccs.${chunk}.bam", overwrite: true
    // publishDir "${params.logs}/ccs", mode: "copy", pattern: "*.log", overwrite: true

    input:
        file(subreads) from raw_subreads_2
        file(index) from raw_index_ch
        each chunk from 1..params.chunks
        val chunks from params.chunks

    output:
        file "*.ccs.${chunk}.bam" into ccs_chunks_ch
        file "ccs.log" into css_log
// need to change ccs.log to ${chunk}.ccs.log
    script:
        """
        ccs \
            --min-rq 0.9 \
            --log-level DEBUG \
            --log-file ccs.log \
            --num-threads ${task.cpus} \
            --chunk ${chunk}/${chunks} \
            ${subreads} \
            ${params.runid}.ccs.${chunk}.bam
        """

}

process ccs_chunk_merging {
    conda "bioconda::pbbam"

    tag "CCS chunk merging"
    // publishDir "${params.ccs}", mode: "copy", pattern: "*.ccs.bam", overwrite: true
    // publishDir "${params.ccs}", mode: "copy", pattern: "*.pbi", overwrite: true

    input:
        file bam_chunks from ccs_chunks_ch.collect()

    output:
        file "*.ccs.bam" into ccs_ch
        file "*.pbi" into index_ch //delete?

    script:
        """
        pbmerge -o ${params.runid}.ccs.bam ${bam_chunks}
        pbindex ${params.runid}.ccs.bam
        """
}

process demux {
    conda "bioconda::lima"

    tag "Demultiplexing samples"
    //publishDir "${params.demux}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.logs}/demux", mode: "copy", pattern: "*.log", overwrite: true
    //publishDir "${params.logs}/demux", mode: "copy", pattern: "*.log", overwrite: true
    input:
        // val sample from sample_name_ch
        file called_ccs from ccs_ch
        file(barcodes) from barcodes_ch
    
    output:
        file "demuxed.*.bam" into demuxed_ch

    script:
        """
        lima \
            --isoseq \
            --dump-clips \
            --peek-guess \
            --num-threads ${task.cpus} \
            --log-level INFO \
            --log-file lima.demux.log \
            --split-bam \
            ${called_ccs} \
            ${barcodes} \
            demuxed.bam
        """
}

process refine {
    conda "bioconda::isoseq3"

    tag "Refining"
    //publishDir "${params.refine}", mode: "copy", pattern: "*.flnc.bam", overwrite: true
    publishDir "${params.logs}/refine", mode: "copy", pattern: "*.log", overwrite: true

    input:
        each file(demuxed) from demuxed_ch.flatten()
        file(barcodes) from refine_barcodes_ch
    
    output:
        file "*.flnc.bam" into refined_ch

    script:
    """
    isoseq3 \
        refine \
        --num-threads ${task.cpus} \
        --log-file ${demuxed.baseName}.log \
        --log-level INFO \
        --verbose \
        ${demuxed} \
        ${barcodes} \
        ${demuxed.baseName}.flnc.bam
    """
}

process cluster {
    conda "bioconda::isoseq3"

    tag "Clustering"
    //publishDir "${params.unpolished}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.logs}/unpolished", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file refined_bam from refined_ch
    
    output:
        file "*.unpolished.bam" into unpolished_bam_ch
        file "*.hq.fasta.gz" into gzipped_hq_ch
        file "*.lq.fasta.gz" into unpolished_lq_ch


    script:
    """
    isoseq3 \
        cluster \
        --use-qvs \
        --num-threads ${task.cpus} \
        --log-file ${refined_bam.baseName}.log \
        --log-level INFO \
        --verbose \
        ${refined_bam} \
        ${refined_bam}.unpolished.bam
    """
}

// Since this takes forever, need to make it optional
// process polish {
//     conda "bioconda::isoseq3"

//     tag "Polishing"
//     //publishDir "${params.polished}", mode: "copy", pattern: "*.bam", overwrite: true
//     publishDir "${params.logs}/polished", mode: "copy", pattern: "*.log", overwrite: true

//     input:
//         // val sample from sample_name_ch
//         file unpolished_bam from unpolished_bam_ch
    
//     output:
//         file "*.polished.bam" into polished_bam_ch
//         file "*.hq.fasta.gz" into polished_hq_ch, gzipped_hq_ch
//         file "*.lq.fasta.gz" into polished_lq_ch

//     script:
//     """
//     isoseq3 \
//         polish \
//         --num-threads ${task.cpus} \
//         --log-file ${unpolished_bam.baseName}.log \
//         --log-level INFO \
//         --verbose \
//         ${unpolished_bam} \
//         ${unpolished_bam.baseName}.polished.bam
//     """
// }

process uncompress {
    conda "bioconda::samtools==1.10"

    tag "Uncompressing"
    //publishDir "${params.transcompressed}", mode: "copy", pattern: "*.fasta.gz", overwrite: true

    input:
        // val sample from sample_name_ch
        file gzipped_fasta from gzipped_hq_ch

    output:
        file "*.fasta" into unpolished_hq_fa_ch, unzipped_hq_fa_ch, filter_fa_hq_ch, unpolished_hq_fa_for_cupcake_ch

    script:
    """
    bgzip --decompress ${gzipped_fasta}
    """
}

process mapping {
    // so I guess gmap 2020.04.08 from bioconda cannot handle compressed fasta?
    conda "bioconda::gmap==2020.04.08"

    tag "Mapping"
    //publishDir "${params.mapped}", mode: "copy", pattern: "*.sam", overwrite: true
    publishDir "${params.logs}/mapped", mode: "copy", pattern: "*.log", overwrite: true

    input:
        file hq_fasta from unpolished_hq_fa_ch

    output:
        file "*.mapped.sam" into mapped_ch

    script:
    """
    gmap \
        -D ${params.reference} \
        -d ${params.species} \
        -f samse \
        -n 0 \
        -t ${task.cpus} \
        --cross-species \
        --max-intronlength-ends 200000 \
        -z sense_force \
        ${hq_fasta} \
        > ${hq_fasta.baseName}.mapped.sam \
        2> ${hq_fasta.baseName}.mapped.log
    """
}

process sort {
    conda "bioconda::samtools==1.10"

    tag "Sorting"
    publishDir "${params.sorted}", mode: "copy", pattern: "*.bam", overwrite: true

    input:
        // val sample from sample_name_ch
        file mapped_sam from mapped_ch

    output:
        file "*.sorted.bam" into sorted_ch, sorted_filtering_bam_ch
    
    script:
    """
    samtools \
        sort \
        -O BAM \
        ${mapped_sam} \
        -o ${mapped_sam.baseName}.sorted.bam
    """
}

process index_sorted {
    conda "bioconda::samtools==1.10"

    tag "Index sorted"

    input:
        file sorted_bam from sorted_ch
    
    output:
        file "*.bai" into sorted_index_ch

    script:
    """
    samtools \
        index \
        -@ ${task.cpus} \
        ${sorted_bam} \
        ${sorted_bam}.bai
    """
}

process filter {
    tag "Filtering"
    //publishDir "${params.filtered}", mode: "copy", pattern: "*.bam", overwrite: true
    
    conda "./filter_sam.yml"

    input:
        // val sample from sample_name_ch
        file hq_fasta from filter_fa_hq_ch
        file sorted from sorted_filtering_bam_ch
        file sorted_index from sorted_index_ch

    output:
        file "${sorted.baseName}_filtered.sam" into filtered_ch

    script:
    """
    filter_sam \
        --fasta ${hq_fasta} \
        --sam ${sorted} \
        --prefix ${sorted.baseName}_
    """
}

process compress {
    conda "bioconda::samtools==1.10"

    tag "Compressing"
    //publishDir "${params.transcompressed}", mode: "copy", pattern: "*.fasta.gz", overwrite: true

    input:
        // val sample from sample_name_ch
        file unzipped_fasta from unzipped_hq_fa_ch

    output:
        file "*.fasta.gz" into bgzipped_hq_ch, bgzipped_hq_ch_2

    script:
    """
    bgzip --index ${unzipped_fasta}
    """
}

process collapse {
    tag "Collapsing"
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.gff", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.unfuzzy", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.group.txt", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.group.txt.unfuzzy", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.rep.fa", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.ignored_ids.txt", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.sam", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.fasta", overwrite: true

    conda "./cdna_cupcake.yml"

    input:
        // val sample from sample_name_ch
        file filtered from filtered_ch
        file fasta from unpolished_hq_fa_for_cupcake_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file "*.collapsed.gff" into collapsed_gff_ch
        // file "*.collapsed.gff.unfuzzy" into collapsed_unfuzzy_gff_ch
        //  file "*.collapsed.rep.fa" into collapsed_fa_ch
    
    script:
    """
    collapse_isoforms_by_sam.py \
        --input ${fasta} \
        --sam ${filtered} \
        --dun-merge-5-shorter \
        --prefix ${fasta.baseName}
    """
}

// fun fact: the perl used by gmst in the SQANTI tool CANNOT tolerate periods in a FASTA sequence name
// guess what EVERY FECKING SEQUENCE NAME HAS TWO OF?
process rename {
    tag "sed name fix"
    conda "./pb.yml"
    
    input:
        file collapsed_gff from collapsed_gff_ch

    output:
        file "*.fixed.gff" into fixed_name_gff_ch

    // two passes because there is a gene_id and transcript_id to fix
    // and while I can make the capture group optional, I cannot seem to
    // figure out how to not add the second period if the transcript number 
    // is missing
    script:
    """
    sed -r 's/PB\.([[:alnum:]]+)\.([[:alnum:]]+)/PB_\1_\2/g' ${collapsed_gff} | \
    sed -r 's/PB\.([[:alnum:]]+)/PB_\1/g' > ${collapsed_gff.baseName}.fixed.gff
    """
}

process sqanti {
    tag "SQANTI3"
    container "milescsmith/sqanti:1.0.0"

    input:
        file fixed_gff from fixed_name_gff_ch
        val chunks from params.chunks

    output:

    script:
    """
    sqanti3_qc.py \
        ${fixed_gff} \
        ${params.annotation} \
        ${params.genome} \
        --cage_peak ${params.cage_peaks} \
        --polyA_motif_list ${params.polyA_list} \
        --cpus {task.cpus} \
        --chunks ${chunks}
    """

}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
