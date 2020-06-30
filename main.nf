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
    .set{ raw_subreads }

Channel
    .value()
    .fromPath( params.barcodes )
    .splitCsv( header:false )
    .map{ row -> row[0].split("-")[0] }
    .unique()
    .set{ sample_name_ch }

extract_bc = { item -> 
    item =~ /bc(\d+)\-\[FR]/
}

process ccs_calling {
    conda "bioconda::pbccs"

    tag "CCS calling"
    publishDir "${params.css}", mode: "copy", pattern: "*.ccs.bam", overwrite: true
    publishDir "${params.css}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        file(subreads) from raw_subreads

    output:
        file "*.ccs.bam" into ccs_ch
        file "ccs.log" into css_log

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

process demux {
    container "quay.io/biocontainers/lima:1.11.0--0"

    cpus 16
    clusterOptions "--mem=256 --parition=highmem -o ~/demux.log"

    tag "Demultiplexing samples"
    publishDir "${params.demux}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.demux}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file called_ccs from ccs_ch
    
    output:
        file "*.bam" into demuxed_ch
        file "*.log" into demux_log

    script:
        """
        lima \
            --isoseq \
            --ccs \
            --dump-clips \
            --peek-guess \
            --num-threads ${task.cpus} \
            --log-level INFO \
            --log-file lima.demux.log \
            --split-bam-named \
            ${called_ccs} \
            ${params.barcodes} \
            demuxed.bam
        """
}

process refine {
    // conda "bioconda::isoseq3"
    container "quay.io/biocontainers/isoseq3:3.3.0--0"

    tag "Refining"
    publishDir "${params.refine}", mode: "copy", pattern: "*.flnc.bam", overwrite: true
    publishDir "${params.refine}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file demuxed_file from demuxed_ch
    
    output:
        file "${refined_bam}" into refined_ch
        file "${refined_log}" into refined_log_ch

    script:
    """
    isoseq3 \
        refine \
        --num-threads ${tasks.cpus} \
        --log-file ${refined_log} \
        --log-level INFO \
        --verbose \
        ${demuxed_file} \
        primers.fasta \
        ${refined_bam}
    """
}

process cluster {
    // conda "bioconda::isoseq3"
    container "quay.io/biocontainers/isoseq3:3.3.0--0"

    tag "Clustering"
    publishDir "${params.unpolished}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.unpolished}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file refined_bam from refined_ch
    
    output:
        file "${unpolished_bam}" into unpolished_bam_ch
        file "${unpolished_log}" into unpolished_log_ch

    script:
    """
    isoseq3 \
        cluster \
        --use-qvs \
        --num-threads ${task.cpu} \
        --log-file ${unpolished_log} \
        --log-level INFO \
        --verbose \
        ${refined_bam} \
        ${unpolished_bam}
    """
}

process polish {
    // conda "bioconda::isoseq3"
    container "quay.io/biocontainers/isoseq3:3.3.0--0"

    tag "Polishing"
    publishDir "${params.polished}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.polished}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file unpolished_bam from unpolished_bam_ch
    
    output:
        file "${polished_bam}" into polished_bam_ch
        file "${polished_log}" into polished_log_ch
        file "*.hq.fasta.gz" into polished_hq_ch, gzipped_hq_ch
        file "*.lq.fasta.gz" into polished_lq_ch

    script:
    """
    isoseq3 \
        polish \
        --num-threads ${task.cpu} \
        --log-file ${polished_log} \
        --log-level INFO \
        --verbose \
        ${unpolished_bam} \
        ${polished_bam}
    """
}

process mapping {
    // conda "bioconda::gmap==2020.04.08"
    container "quay.io/biocontainers/gmap:2020.06.01--pl526h2f06484_1"

    tag "Mapping"
    publishDir "${params.mapped}", mode: "copy", pattern: "*.sam", overwrite: true
    publishDir "${params.mapped}", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file polished_hq_fasta from polished_hq_ch

    output:
        file "${mapped_sam}" into mapped_ch
        file "${mapping_log}" into mapping_log_ch

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
        ${polished_hq_fasta} \
        > ${mapped_sam} \
        2> ${mapping_log}
    """
}

process sort {
    // conda "bioconda::samtools==1.10"
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    tag "Sorting"
    publishDir "${params.sorted}", mode: "copy", pattern: "*.bam", overwrite: true

    input:
        // val sample from sample_name_ch
        file mapped_sam from mapped_ch

    output:
        file "${sorted_bam}" into sorted_ch
    
    script:
    """
    samtools \
        sort \
        -O BAM \
        ${mapped_sam} \
        -o ${sorted_bam}
    """
}

process transcompress {
    // conda "bioconda::samtools==1.10"
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    tag "Transcompressing"
    publishDir "${params.transcompressed}", mode: "copy", pattern: "*.fasta.gz", overwrite: true

    input:
        // val sample from sample_name_ch
        file "${sample}.fasta.gz" from gzipped_hq_ch

    output:
        file "*.fasta.gz" into bgzipped_hq_ch, bgzipped_hq_ch_2

    script:
    """
    bgzip \
        --decompress ${sample}.fasta.gz && \
    bgzip \
        --index ${sample}.fasta
    """
}

process filter {
    container "registry.gitlab.com/milothepsychic/filter_sam"
    
    tag "Filtering"
    publishDir "${params.filtered}", mode: "copy", pattern: "*.bam", overwrite: true

    input:
        // val sample from sample_name_ch
        file "${sample}.fasta.gz" from bgzipped_hq_ch
        file sorted from sorted_ch

    output:
        file "${sample}_filtered.sam" into filtered_ch

    script:
    """
    filter_sam \
        --fasta ${sample}.fasta.gz \
        --sam ${sorted} \
        --prefix ${sample}_
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

    container "milescsmith/cdna_cupcake:./cdna_cupcake.yml"

    input:
        // val sample from sample_name_ch
        file filtered from filtered_ch
        file fasta from unpolished_hq_fa_for_cupcake_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file "*.collapsed.gff" into collapsed_gff_ch
        // file "*.collapsed.gff.unfuzzy" into collapsed_unfuzzy_gff_ch
        file "*.collapsed.rep.fa" into collapsed_fa_ch
    
    script:
    """
    collapse_isoforms_by_sam.py \
        --input ${fasta} \
        --sam ${filtered} \
        --dun-merge-5-shorter \
        --prefix ${fasta.baseName}
    """
}

// perl - which is unfortunately still key to GeneMark S-T use
// by pygmst - seems to have a problem with either double dashes
// or long file names.  So we are going to rename everything
process rename {
    tag "sed name fix"
    container "milescsmith/rename:0.20-7"
    
    input:
        file collapsed_fa from collapsed_fa_ch

    output:
        file "*.fa" into fixed_name_fa_ch

    // two passes because there is a gene_id and transcript_id to fix
    // and while I can make the capture group optional, I cannot seem to
    // figure out how to not add the second period if the transcript number 
    // is missing
    script:
    """
    rename 's/_5p\-\-bc[0-9]{4}_3p//' ${collapsed_fa} | rename 's/^demuxed\.//'
    """
}

process sqanti {
    tag "SQANTI3"
    container "milescsmith/sqanti:1.3.11"

    publishDir "${params.sqanti}", mode: "copy", pattern: "*.pdf", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep.params.txt", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_classification.txt", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.faa", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.fasta", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.fasta.fai", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.genePred", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.gtf", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.gtf.cds.gff", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.gtf.tmp", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected.sam", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_corrected_indels.txt", overwrite: true
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.rep_junctions.txt", overwrite: true

    input:
        file fixed_name_fa from fixed_name_fa_ch

    output:
        // *.rep.params.txt
        // *.rep_classification.txt
        // *.rep_corrected.faa
        // *.rep_corrected.fasta
        // *.rep_corrected.fasta.fai
        // *.rep_corrected.genePred
        // *.rep_corrected.gtf
        // *.rep_corrected.gtf.cds.gff
        // *.rep_corrected.gtf.tmp
        // *.rep_corrected.sam
        // *.rep_corrected_indels.txt
        // *.rep_junctions.txt
        file "*.rep_sqanti_report.pdf" into sq_report_ch

    script:
    """
    sqanti3_qc \
        /s/guth-aci/isoseq/11_collapsed/*.rep.fa \
        ${params.annotation} \
        ${params.genome} \
        --cage_peak ${params.cage_peaks} \
        --polyA_motif_list ${polyA} \
        --cpus ${task.cpus}
    """
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
