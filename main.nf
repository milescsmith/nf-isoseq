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
// Date: 2021/05/05
// Version: 0.2.0

// Note: I've borrowed a bit of code from nf-core/rnaseq

// File locations and program parameters
nextflow.enable.dsl=2

checkPathParamList = [
    params.annotation,
    params.genome,
    params.cage_peaks,
    params.polyA_list,
    params.barcode,
]

alignerList = ["minimap2", "desalt2", "gmap"]

for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
        } 
}

if (!alignerList.contains(params.aligner)) {
    exit 1, "That is not a valid aligner choice.  Please choose from ${alignerList.join(', ')}"
}

def helpMessage() {
    log.info nfcoreHeader()
    log.info """
    Usage:

      The typical command for running the pipeline is as follows:
 
      nextflow run milescsmith/nf-isoseq \\
        --project /path/to/project \\
        --genome /path/to/reference/genome.fa \\
        --annotation /path/to/reference/annotation.gff \\
        --cage_preaks /path/to/reference/cage_peaks.bed \\
        --polyA_list /path/to/reference/poly_a_list.txt \\
        --barcodes /path/to/barcodes \\
        --profile slurm

    Mandatory arguments:
      --project             Directory to use as a base where raw files are 
                            located and results and logs will be written
      --profile             Configuration profile to use for processing.
                            Available: slurm, gcp, standard
      --species
      --reference
      --annoation
      --genome
      --cage_peaks
      --polyA_list
      --barcodes

    Processing options:
      --chunks              Number of chunks to divide the initial BAM file into
                             processing by ccs

    Results locations:      If not specified, these will be created relative to
                            the `project` argument
                            Note: if undefined, the pipeline expects the raw BAM
                            file to be located in `/path/to/project/01_raw`
      --logs                Default: project/logs
      --raw                 Default: project/01_raw
      --ccs                 Default: project/02_ccs
      --demux               Default: project/03_demux
      --refined             Default: project/04_refined
      --unpolished          Default: project/05_unpolished
      --polished            Default: project/06_polished
      --mapped              Default: project/07_mapped
      --sorted              Default: project/08_sorted
      --transcompressed     Default: project/09_transcompress
      --filtered            Default: project/10_filtered
      --collapsed           Default: project/11_collapsed
      --sqanti              Default: project/12_sqanti
    
    Reference locations:

    
    Other:
      --runid               Prefix to use in naming certain files
      --help                Show this message
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


params.input = "${params.raw}/*.subreads.bam"

Channel
    .fromPath( params.input , checkIfExists: true )
    .into{ raw_subreads_1; raw_subreads_2 }

// Channel
//     .value()
//     .fromPath( params.barcodes )
//     .splitCsv( header:false )
//     .map{ row -> row[0].split("-")[0] }
//     .unique()
//     .set{ sample_name_ch }

Channel
    .fromPath( params.barcodes )
    .into{ barcodes_ch; refine_barcodes_ch }

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

// need to figure out how to rename the files made here
// something in GMST does NOT like dashes in the file name
process demux {
    conda "bioconda::lima"
    // container "quay.io/biocontainers/lima:1.11.0--0"

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
    // container "quay.io/biocontainers/isoseq3:3.3.0--0"

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
    // container "quay.io/biocontainers/isoseq3:3.3.0--0"

    tag "Clustering"
    //publishDir "${params.unpolished}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.logs}/unpolished", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file refined_bam from refined_ch
    
    output:
        file "*.unpolished.bam" into unpolished_bam_ch
        // file "*.hq.fasta.gz" into gzipped_hq_ch
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

//Since this takes forever, need to make it optional
process polish {
    conda "bioconda::isoseq3"
// container "quay.io/biocontainers/isoseq3:3.3.0--0"

    tag "Polishing"
    //publishDir "${params.polished}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir "${params.logs}/polished", mode: "copy", pattern: "*.log", overwrite: true

    input:
        // val sample from sample_name_ch
        file unpolished_bam from unpolished_bam_ch
    
    output:
        file "*.polished.bam" into polished_bam_ch
        file "*.hq.fasta.gz" into polished_hq_ch, gzipped_hq_ch
        file "*.lq.fasta.gz" into polished_lq_ch

    script:
    """
    isoseq3 \
        polish \
        --num-threads ${task.cpus} \
        --log-file ${unpolished_bam.baseName}.log \
        --log-level INFO \
        --verbose \
        ${unpolished_bam} \
        ${unpolished_bam.baseName}.polished.bam
    """
}

process uncompress {
    conda "bioconda::samtools==1.10"
    // container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    tag "Uncompressing"
    //publishDir "${params.transcompressed}", mode: "copy", pattern: "*.fasta.gz", overwrite: true

    input:
        // val sample from sample_name_ch
        file gzipped_fasta from gzipped_hq_ch

    output:
        file "*.fasta" into unpolished_hq_fa_ch, unpolished_hq_fa_for_cupcake_ch //, unzipped_hq_fa_ch, filter_fa_hq_ch

    script:
    """
    bgzip --decompress ${gzipped_fasta}
    """
}

process mapping {
    // so I guess gmap 2020.04.08 from bioconda cannot handle compressed fasta?
    conda "bioconda::gmap==2020.06.01"
    // container "quay.io/biocontainers/gmap:2020.06.01--pl526h2f06484_1"

    tag "Mapping"
    // publishDir "${params.mapped}", mode: "copy", pattern: "*.sam", overwrite: true
    publishDir "${params.logs}/mapped", mode: "copy", pattern: "*.log", overwrite: true

    input:
        file hq_fasta from unpolished_hq_fa_ch

    output:
        tuple file("*.mapped.sam"), file(hq_fasta) into mapped_ch

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
    // container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    tag "Sorting"
    publishDir "${params.sorted}", mode: "copy", pattern: "*.bam", overwrite: true

    input:
        // val sample from sample_name_ch
        tuple file(mapped_sam), file(hq_fasta) from mapped_ch

    output:
        tuple file("*.sorted.sam"), file(hq_fasta) into sorted_for_collapsing_bam_ch, sorting_ch //sorted_filtering_bam_ch
    
    script:
    """
    samtools \
        sort \
        -O SAM \
        ${mapped_sam} \
        -o ${mapped_sam.baseName}.sorted.sam
    """
}

// process index_sorted {
//     conda "bioconda::samtools==1.10"
//     // container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

//     tag "Index sorted"

//     input:
//         tuple file(sorted_sam), file(hq_fasta) from sorting_ch
    
//     output:
//         file "*.bai" into sorted_index_ch

//     script:
//     """
//     samtools \
//         index \
//         -@ ${task.cpus} \
//         ${sorted_sam} \
//         ${sorted_sam}.bai
//     """
// }

// process filter {
//     tag "Filtering"
//     //publishDir "${params.filtered}", mode: "copy", pattern: "*.bam", overwrite: true
    
//     conda "./environments/filter_sam.yml"
//     // container "registry.gitlab.com/milothepsychic/filter_sam"

//     input:
//         // val sample from sample_name_ch
//         file hq_fasta from filter_fa_hq_ch
//         file sorted from sorted_filtering_bam_ch
//         file sorted_index from sorted_index_ch

//     output:
//         file "${sorted.baseName}_filtered.sam" into filtered_ch

//     script:
//     """
//     filter_sam \
//         --fasta ${hq_fasta} \
//         --sam ${sorted} \
//         --prefix ${sorted.baseName}_
//     """
// }

// process compress {
//     conda "bioconda::samtools==1.10"
//     // container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

//     tag "Compressing"
//     //publishDir "${params.transcompressed}", mode: "copy", pattern: "*.fasta.gz", overwrite: true

//     input:
//         // val sample from sample_name_ch
//         file unzipped_fasta from unzipped_hq_fa_ch

//     output:
//         file "*.fasta.gz" into bgzipped_hq_ch, bgzipped_hq_ch_2

//     script:
//     """
//     bgzip --index ${unzipped_fasta}
//     """
// }

process collapse {
    tag "cDNA_Cupcake Collapse"
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.gff", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.unfuzzy", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.group.txt", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.group.txt.unfuzzy", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.rep.fa", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.ignored_ids.txt", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.sam", overwrite: true
    publishDir "${params.collapsed}", mode: "copy", pattern: "*.fasta", overwrite: true

    conda "./environments/cdna_cupcake.yml"
    // container "milescsmith/cdna_cupcake:12.2.8"

    input:
        // val sample from sample_name_ch
        // file filtered from filtered_ch
        tuple file(sorted_sam), file(hq_fasta) from sorted_for_collapsing_bam_ch
        // file fasta from unpolished_hq_fa_for_cupcake_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file "*.collapsed.gff" into collapsed_gff_ch
        // file "*.collapsed.gff.unfuzzy" into collapsed_unfuzzy_gff_ch
        file "*.collapsed.rep.fa" into collapsed_fa_ch
    
    script:
    """
    collapse_isoforms_by_sam \
        --input ${hq_fasta} \
        --sam ${sorted_sam} \
        --dun-merge-5-shorter \
        --prefix ${hq_fasta.baseName}
    """
}

// GeneMark S-T which is still foundational to pygmst - seems
// to have a problem with either double dashes
// or long file names.  So we are going to rename everything
process rename {
    tag "name fix GeneMark"
    container "milescsmith/rename:0.20-7"
    
    input:
        file collapsed_fa from collapsed_fa_ch

    output:
        file "*.fa" into fixed_name_fa_ch

    script:
    """
    rename 's/_5p\\-\\-bc[0-9]{4}_3p//' ${collapsed_fa} | rename 's/^demuxed\\.//'
    """
}

process sqanti {
    tag "SQANTI3"
    container "milescsmith/sqanti3:1.4.0"
    errorStrategy 'ignore' // sometimes, there's just this one file...
    // conda "./environments/sqanti3.yml"

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
    publishDir "${params.sqanti}", mode: "copy", pattern: "*.unpolished.hq.collapsed.rep.genePred", overwrite: true

    input:
        file fixed_name_fa from fixed_name_fa_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file "*.rep.params.txt" into sqanti_ch_1
        file "*.rep.renamed.fasta" into sqanti_ch_2
        file "*.rep_classification.txt" into sqanti_ch_3
        file "*.rep_corrected.faa" into sqanti_ch_4
        file "*.rep_corrected.fasta" into sqanti_ch_5
        file "*.rep_corrected.fasta.fai" into sqanti_ch_6
        file "*.rep_corrected.genePred" into sqanti_ch_7
        file "*.rep_corrected.gtf" into sqanti_ch_8
        file "*.rep_corrected.gtf.cds.gff" into sqanti_ch_9
        file "*.rep_corrected.gtf.tmp" into sqanti_ch_10
        file "*.rep_corrected.sam" into sqanti_ch_11
        file "*.rep_corrected_indels.txt" into sqanti_ch_12
        file "*.rep_junctions.txt" into sqanti_ch_13
        file "*.unpolished.hq.collapsed.rep.genePred" into sqanti_ch_14
        file "*.rep_sqanti_report.pdf" into sq_report_ch

    script:
    """
    sqanti3_qc \
        ${fixed_name_fa} \
        ${params.annotation} \
        ${params.genome} \
        --cage_peak ${params.cage_peaks} \
        --polyA_motif_list ${params.polyA_list} \
        --cpus ${task.cpus}
    """
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnaseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
