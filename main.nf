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

// TODO: convert to DSL2
// TODO: check parameters on run to make sure they are valid
// TODO: check all files present when starting
// TODO: need functionality to download/build indices, grab cage_peaks and polyA lists
// TODO: enable more options

// Nextflow pipeline for processing PacBio IsoSeq runs
// Author: Miles Smith <miles-smith@omrf.org>
// Date: 2020/05/26
// Version: 0.1.0

// File locations and program parameters

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
    .fromPath( params.input, checkIfExists: true )
    .into{ raw_subreads_1; raw_subreads_2 }

Channel
    .fromPath( params.barcodes, checkIfExists: true )
    .set{ barcodes_ch }

// Channel
//     .fromPath( params.barcodes, checkIfExists: true )
//     .set{ refine_barcodes_ch }

process ccs_indexing {
    tag "CCS indexing"
    // publishDir "${params.raw}",
        // mode: "copy",
        // pattern : "*.bam.pbi",
        // overwrite: false
    
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
    
    tag "CCS calling"
    
    publishDir "${params.logs}/ccs", mode: "copy", pattern: "*.txt", overwrite: true
    publishDir "${params.logs}/ccs", mode: "copy", pattern: "*.json.gz", overwrite: true

    input:
        file(subreads) from raw_subreads_2
        file(index)    from raw_index_ch
        each(chunk)    from 1..params.chunks
        val(chunks)    from params.chunks

    output:
        file("*.bam")            into ccs_chunks_ch
        file("*.ccs_report.txt") into ccs_log_ch
        file("*.json.gz")        into ccs_metrics_ch

    script:
        if (params.runid != "")
            """
            ccs \
                --min-rq 0.9 \
                --log-level ${params.calling_log_level} \
                --log-file ccs.${chunk}.log \
                --num-threads ${task.cpus} \
                --chunk ${chunk}/${chunks} \
                ${subreads} \
                ${params.runid}.ccs.${chunk}.bam
            """
        else
            """
            ccs \
                --min-rq 0.9 \
                --log-level ${params.calling_log_level} \
                --log-file ccs.${chunk}.log \
                --num-threads ${task.cpus} \
                --chunk ${chunk}/${chunks} \
                ${subreads} \
                ccs.${chunk}.bam
            """

}

process ccs_chunk_merging {

    tag "CCS chunk merging"

    input:
        file(bam_chunks) from ccs_chunks_ch.collect()

    output:
        file("*.ccs.bam") into ccs_ch
        file("*.pbi") into index_ch //delete?

    script:
        if (params.runid != "")
            """
            pbmerge -o ${params.runid}.ccs.bam ${bam_chunks}
            pbindex ${params.runid}.ccs.bam
            """
        else
            """
            pbmerge -o merged.ccs.bam ${bam_chunks}
            pbindex merged.ccs.bam
            """
}

process demux { 

    tag "Demultiplexing samples"

    publishDir path: "${params.logs}/demux", mode: "copy", pattern: "*.log", overwrite: true
    
    input:
        // val sample from demux_sample_name_ch
        file(called_ccs) from ccs_ch
        file(barcodes)   from barcodes_ch

    output:
        file("demuxed.*.bam")  into demuxed_ch
        file("lima.demux.log") into demuxed_log_ch

    script:
        if (params.biosample != "")
            """
            lima \
                --isoseq \
                --peek-guess \
                --num-threads ${task.cpus} \
                --log-level ${params.demux_log_level} \
                --log-file lima.demux.log \
                --split-bam-named \
                --biosample-csv ${params.biosample} \
                ${called_ccs} \
                ${barcodes} \
                demuxed.bam
            """
        else
            """
            lima \
                --isoseq \
                --peek-guess \
                --num-threads ${task.cpus} \
                --log-level ${params.demux_log_level} \
                --log-file lima.demux.log \
                --split-bam \
                --split-bam-named \
                ${called_ccs} \
                ${barcodes} \
                demuxed.bam
            """

}


process refine {
    
    tag "Refining"

    publishDir path: "${params.refined}/${sample_name}", mode: "copy", pattern: "*.bam", overwrite: true, saveAs: { filename -> "demuxed.flnc.bam" }
    publishDir path: "${params.refined}/${sample_name}", mode: "copy", pattern: "*.bam.pbi", overwrite: true, saveAs: { filename -> "demuxed.flnc.bam.pbi" }
    publishDir path: "${params.refined}/${sample_name}", mode: "copy", pattern: "*.xml", overwrite: true, saveAs: { filename -> "demuxed.flnc.consensusreadset.xml" }
    publishDir path: "${params.refined}/${sample_name}", mode: "copy", pattern: "*.json", overwrite: true, saveAs: { filename -> "demuxed.flnc.filter_summary.json" }
    publishDir path: "${params.refined}/${sample_name}", mode: "copy", pattern: "*.csv", overwrite: true, saveAs: { filename -> "demuxed.post_refine_report.csv" }
    publishDir path: "${params.logs}/refine", mode: "copy", pattern: "*.log", overwrite: true

    input:
        file(demuxed) from demuxed_ch.flatten()
        // file(barcodes) from refine_barcodes_ch
    
    output:
        tuple val(sample_name), file("*.bam") into refined_ch
        file("*.csv") into refined_report_ch
        file("*.log") into refined_log_ch
        file("*.xml") into refined_consensusreadset_ch
        file("*.pbi") into refined_index_ch
        file("*.json") into refined_filter_ch

    script:
        sample_name = (demuxed.baseName =~ /demuxed\.([\s\S]*)_5p--[\s\S]*_3p/)[0][1]
        """
        isoseq3 \
            refine \
            --num-threads ${task.cpus} \
            --log-file refining.${sample_name}.log \
            --log-level ${params.refining_log_level} \
            --verbose \
            ${demuxed} \
            ${params.barcodes} \
            demuxed.refined.${sample_name}.flnc.bam
        """
}

process cluster {

    tag "Clustering"

    publishDir path: "${params.unpolished}/${sample_name}", mode: "copy", pattern: "*.hq.bam", overwrite: true, saveAs: { filename -> "flnc.unpolished.hq.bam" }
    publishDir path: "${params.unpolished}/${sample_name}", mode: "copy", pattern: "*.lq.bam", overwrite: true, saveAs: { filename -> "flnc.unpolished.lq.bam" }
    publishDir path: "${params.unpolished}/${sample_name}", mode: "copy", pattern: "*.fasta.gz", overwrite: true, saveAs: { filename -> "flnc.unpolished.fasta.gz" }
    publishDir path: "${params.unpolished}/${sample_name}", mode: "copy", pattern: "*.csv", overwrite: true, saveAs: { filename -> "flnc.unpolished.cluster_report.csv" }
    publishDir path: "${params.unpolished}/${sample_name}", mode: "copy", pattern: "*.cluster", overwrite: true, saveAs: { filename -> "flnc.unpolished.cluster" }
    publishDir path: "${params.logs}/unpolished", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "clustering.${sample_name}.log" }

    input:
        tuple val(sample_name), file(refined_bam) from refined_ch
    
    output:
        tuple val(sample_name), file("*.unpolished.bam") into unpolished_bam_ch
        file("*.hq.fasta.gz") into gzunpolished_hq_ch
        file("*.lq.fasta.gz") into gzunpolished_lq_ch
        file("*.csv") into cluster_report_ch
        file("*.cluster") into clusters_ch
        file("*.log") into cluster_log_ch


    script:
        """
        isoseq3 \
            cluster \
            --use-qvs \
            --num-threads ${task.cpus} \
            --log-file ${refined_bam.baseName}.log \
            --log-level ${params.cluster_log_level} \
            --verbose \
            ${refined_bam} \
            ${refined_bam}.unpolished.bam
        """

}

// Since this takes forever, need to make it optional or find a way to speed it up
process polish {
    
    tag "Polishing"

    publishDir path: "${params.polished}/${sample_name}", mode: "copy", pattern: "*.bam", overwrite: true, saveAs: { filename -> "demuxed.polished.bam"}
    publishDir path: "${params.polished}/${sample_name}", mode: "copy", pattern: "*.hq.fasta.gz", overwrite: true, saveAs: {filename -> "demuxed.hq.fasta.gz"}
    publishDir path: "${params.polished}/${sample_name}", mode: "copy", pattern: "*.lq.fasta.gz", overwrite: true, saveAs: {filename -> "demuxed.lq.fasta.gz"}
    publishDir path: "${params.logs}/polished", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "polish.${sample_name}.log" }

    input:
        tuple val(sample_name), file(unpolished_bam) from unpolished_bam_ch
        file(hq_fasta) from gzunpolished_hq_ch
        file(lq_fasta) from gzunpolished_lq_ch
    
    output:
        tuple val(sample_name), file("*.polished.bam") into polished_bam_ch
        tuple val(sample_name), file("*.hq.fasta.gz") into polished_hq_ch
        file("*.lq.fasta.gz") into polished_lq_ch
        file("*.log")         into polished_log_ch

    script:
        if (params.polish)
            """
            isoseq3 \
                polish \
                --num-threads ${task.cpus} \
                --log-file ${unpolished_bam.baseName}.log \
                --log-level ${params.polish_log_level} \
                --verbose \
                ${unpolished_bam} \
                ${unpolished_bam.baseName}.polished.bam
            """
        else
            """
            mv ${unpolished_bam} demuxed.polished.bam
            mv ${hq_fasta} demuxed.hq.fasta.gz
            mv ${lq_fasta} demuxed.lq.fasta.gz
            touch ${unpolished_bam.baseName}.log
            """

}

// process decompress {
//     conda "bioconda::samtools==1.10"
//     // container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

//     tag "Decompressing"
//     //publishDir "${params.transcompressed}",
        // mode: "copy",
        // pattern: "*.fasta.gz",
        // overwrite: true

//     input:
//         // val sample from sample_name_ch
//         file gzipped_fasta from gzipped_hq_ch

//     output:
//         file "*.fasta" into polished_hq_fa_ch, unpolished_hq_fa_for_cupcake_ch //, unzipped_hq_fa_ch, filter_fa_hq_ch

//     script:
//     """
//     bgzip --decompress ${gzipped_fasta}
//     """
// }

// Rewrite to make it easy to switch aligners.  Maybe only easy with DSL2
// process mapping {
//     // so I guess gmap 2020.04.08 from bioconda cannot handle compressed fasta?
//     conda "bioconda::gmap==2020.06.01"
//     // container "quay.io/biocontainers/gmap:2020.06.01--pl526h2f06484_1"

//     tag "Mapping"
//     //publishDir "${params.mapped}",
        // mode: "copy",
        // pattern: "*.sam",
        // overwrite: true
//     publishDir "${params.logs}/mapped",
        // mode: "copy",
        // pattern: "*.log",
        // overwrite: true

//     input:
//         file hq_fasta from unpolished_hq_fa_ch

//     output:
//         tuple file("*.mapped.sam"), file(hq_fasta) into mapped_ch

//     script:
//     """
//     gmap \
//         -D ${params.reference} \
//         -d ${params.species} \
//         -f samse \
//         -n 0 \
//         -t ${task.cpus} \
//         --cross-species \
//         --max-intronlength-ends 200000 \
//         -z sense_force \
//         ${hq_fasta} \
//         > ${hq_fasta.baseName}.mapped.sam \
//         2> ${hq_fasta.baseName}.mapped.log
//     """
// }

process mapping {

    tag "Mapping"

    publishDir path: "${params.mapped}/${sample_name}", mode: "copy", pattern: "*.sam", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped.sam" : "hq_unpolished_mapped.sam" }
    publishDir path: "${params.logs}/mapped", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "mapping.${sample_name}.log" }

    input:
        tuple val(sample_name), file(polished_hq_fasta) from polished_hq_ch

    output:
        tuple val(sample_name), file("*.mapped.sam"), file(polished_hq_fasta) into mapped_ch
        file("*.mapped.log") into mapped_log_ch

    script:
        """
        minimap2 \
            -H \
            -t ${task.cpus} \
            -ax splice:hq \
            -uf \
            --secondary=no \
            -O6,24 \
            -B4 \
            ${params.minimap_index} \
            ${polished_hq_fasta} \
            > ${polished_hq_fasta.baseName}.mapped.sam \
            2> ${polished_hq_fasta.baseName}.mapped.log
        """
}

process sort {

    tag "Sorting"

    publishDir path: "${params.sorted}/${sample_name}", mode: "copy", pattern: "*.bam", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted.bam" : "hq_unpolished_mapped_sorted.bam" }
    // publishDir path: "${params.sorted}/${sample_name}", mode: "copy", pattern: "*.sam", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted.${filename.Extension}" : "hq_unpolished_mapped_sorted.${filename.Extension}" }
    publishDir path: "${params.sorted}/${sample_name}", mode: "copy", pattern: "*.fasta", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted.fasta" : "hq_unpolished_mapped_sorted.fasta}" }

    input:
        // val sample from sort_sample_name_ch
        tuple val(sample_name), file(mapped_sam), file(hq_fasta) from mapped_ch

    output:
        tuple val(sample_name), file("*.sorted.sam"), file(hq_fasta) into sorted_for_collapsing_bam_ch, sorting_ch //sorted_filtering_bam_ch
    
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
//     //publishDir "${params.filtered}",
        // mode: "copy",
        // pattern: "*.bam",
        // overwrite: true
    
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
//     //publishDir "${params.transcompressed}",
        // mode: "copy",
        // pattern: "*.fasta.gz",
        // overwrite: true

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
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.gff", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.gff" : "hq_unpolished_mapped_sorted_collapsed.gff" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.unfuzzy", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.unfuzzy" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.group.txt" : "hq_unpolished_mapped_sorted_collapsed.group.txt" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt.unfuzzy", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.group.txt.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.group.txt.unfuzzy" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.rep.fa", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.rep.fa" : "hq_unpolished_mapped_sorted_collapsed.rep.fa" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.ignored_ids.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.ignored_ids.txt" : "hq_unpolished_mapped_sorted_collapsed.ignored_ids.txt" }

    input:
        tuple val(sample_name), file(sorted_sam), file(hq_fasta) from sorted_for_collapsing_bam_ch
        // file filtered from filtered_ch
        // file fasta from unpolished_hq_fa_for_cupcake_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file("*.collapsed.gff")     into collapsed_gff_ch
        tuple val(sample_name), file("*.collapsed.rep.fa")  into collapsed_fa_ch
        file("*.unfuzzy")           into collapsed_unfuzzy_ch
        tuple val(sample_name), file("*.group.txt")         into collapsed_group_ch
        file("*.group.txt.unfuzzy") into collapsed_group_unfuzzy_ch
        file("*.ignored_ids.txt")   into collapsed_ignored_ids_ch
    
    script:
    """
    collapse_isoforms_by_sam \
        --input ${hq_fasta} \
        --sam ${sorted_sam} \
        --dun-merge-5-shorter \
        --prefix ${hq_fasta.baseName}
    """
}


process sqanti {
    
    tag "SQANTI3"
    // errorStrategy 'ignore' // sometimes, there's just this one file...

    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.params.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep.params.txt": "hq_unpolished_mapped_sorted_collapsed_rep.params.txt" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep.renamed.fasta", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_renamed.fasta": "hq_unpolished_mapped_sorted_collapsed_rep_renamed.fasta" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_classification.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_classification.txt": "hq_unpolished_mapped_sorted_collapsed_rep_classification.txt" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.faa", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.faa": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.faa" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.fnn", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fnn": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fnn" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.fasta", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fasta": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fasta" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.fasta.fai", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fasta.fai": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fasta.fai" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.genePred", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.genePred": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.genePred" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.gtf", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.gtf": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.gtf" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.cds.gff", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.cds.gff": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.cds.gff" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected.sam", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected.sam": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.sam" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_corrected_indels.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_corrected_indels.txt": "hq_unpolished_mapped_sorted_collapsed_rep_corrected_indels.txt" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep_junctions.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep_junctions.txt": "hq_unpolished_mapped_sorted_collapsed_rep_junctions.txt" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "*.rep.genePred", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep.genePred": "hq_unpolished_mapped_sorted_collapsed_rep.genePred" }
    publishDir path: "${params.sqanti}/${sample_name}", mode: "copy", pattern: "report.pdf", overwrite: true

    input:
        tuple val(sample_name), file(fixed_name_fa) from collapsed_fa_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file "*.params.txt" into sqanti_ch_1
        file "*.rep.renamed.fasta" into sqanti_ch_2
        file "*.rep_classification.txt" into sqanti_ch_3
        file "*.rep_corrected.faa" into sqanti_ch_4
        file "*.rep_corrected.fnn" into sqanti_ch_5
        file "*.rep_corrected.fasta" into sqanti_ch_6
        file "*.rep_corrected.fasta.fai" into sqanti_ch_7
        file "*.rep_corrected.genePred" into sqanti_ch_8
        file "*.rep_corrected.gtf" into sqanti_ch_9
        file "*.rep_corrected.cds.gff" into sqanti_ch_10
        file "*.rep_corrected.sam" into sqanti_ch_11
        file "*.rep_corrected_indels.txt" into sqanti_ch_12
        file "*.rep_junctions.txt" into sqanti_ch_13
        file "*.rep.genePred" into sqanti_ch_14
        file "report.pdf" into sq_report_ch

    script:
    """
    sqanti3_qc \
        ${fixed_name_fa} \
        ${params.annotation} \
        ${params.genome} \
        --cage_peak ${params.cage_peaks} \
        --polyA_motif_list ${params.polyA_list} \
        --cpus ${task.cpus} \
        -vvv
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
    c_reset  = params.monochrome_logs ? '': "\033[0m";
    c_dim    = params.monochrome_logs ? '': "\033[2m";
    c_black  = params.monochrome_logs ? '': "\033[0;30m";
    c_green  = params.monochrome_logs ? '': "\033[0;32m";
    c_yellow = params.monochrome_logs ? '': "\033[0;33m";
    c_blue   = params.monochrome_logs ? '': "\033[0;34m";
    c_purple = params.monochrome_logs ? '': "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '': "\033[0;36m";
    c_white  = params.monochrome_logs ? '': "\033[0;37m";

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
