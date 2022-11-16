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

    Processing options:
      --chunks              Number of chunks to divide the initial BAM file into
                            processing by ccs

    Results locations:      If not specified, these will be created relative to
                            the `project` argument
                            Note: if undefined, the pipeline expects the raw BAM
                            file to be located in `/path/to/project/01_raw`
      --logs                Default: project/logs
      --raw                 Default: project/01_raw
      --input               
    
    Reference locations:
      --species             Default: homo_sapiens
      --references          Default: /Volumes/guth_aci_informatics/references
      --annotation          Default: references/genomic/${params.species}/sequences/gencode_v38/gencode.v38.primary_assembly.annotation.gtf
      --genome              Default: references/genomic/${params.species}/sequences/gencode_v38/GRCh38.p13.genome.fa
      --cage_peaks          Default: references/miscellaneous/hg38.cage_peak_phase1and2combined_coord.bed
      --polyA_list          Default: references/miscellaneous/human.polyA.list.txt
      --barcodes            Default: project/barcodes.fa
      --gmap_index          Default: references/genomic/${params.species}/indices/gmap/gencode_v32
      --minimap_index       Default: references/genomic/${params.species}/indices/minimap2/gencode_v38/gencode_v38.mmi
      --tappas_annotation   Default: references/genomic/$species/sequences/Homo_sapiens_GRCh38_Ensembl_86_tappAS.gff3

    
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

params.input = "${params.flncs}/flnc-*.bam"

Channel
    .fromPath( params.input, checkIfExists: true )
    .into{ flncs_ch }


process mapping {

    tag "Mapping"

    publishDir path: "${params.mapped}/${sample_name}", mode: "copy", pattern: "*.bam", overwrite: true, saveAs: { filename -> "mapped.${sample_name}.bam" }
    publishDir path: "${params.logs}/mapped",           mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "mapping.${sample_name}.log" }

    input:
        file(input_reads) from flncs_ch

    output:
        tuple val(sample_name), file("*.sam"), file("*.hq.fasta.gz") into mapped_ch
        file("*.alignment.log")                                      into mapped_log_ch

    script:
        sample_name = (input_reads.baseName =~ /flnc\-([\d]*)/)[0][1]
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
            ${input_reads} \
            > sample.${sample_name}.mapped.sam \
            2> sample.${sample_name}.alignment.log
        """

}


process sort {

    tag "Sorting"

    input:
        // val sample from sort_sample_name_ch
        tuple val(sample_name), file(mapped_sam), file(hq_fasta)     from mapped_ch

    output:
        tuple val(sample_name), file("*.sorted.sam")                 into sorted_for_compressing_bam_ch

    script:
        """
        samtools \
            sort \
            -O SAM \
            --reference ${params.ref_seq} \
            ${mapped_sam} \
            -o sample.${sample_name}.sorted.sam
        """

}


process sam_to_bam {

    tag "Compressing"
    publishDir "${params.sorted}", mode: "copy", pattern: "*.bam", overwrite: true

    input:
        tuple val(sample_name), file(sorted_sam) from sorted_for_compressing_bam_ch

    output:
        tuple val(sample_name), file("*.bam")    into sorted_bam_ch

    script:
        """
        samtools \
            view \
            -b \
            -@ ${task.cpus} \
            --reference ${params.ref_seq} \
            ${sorted_sam} \
            -o sample.${sample_name}.mapped.sorted.bam
        """

}


process index_bam {

    tag "Indexing mapped"
    publishDir "${params.sorted}", mode: "copy", pattern: "*.bai", overwrite: true

    input:
        tuple val(sample_name), file(sorted_sam) from sorted_bam_ch

    output:
        tuple val(sample_name), file("*.bai")    into indexed_bam_ch

    script:
        """
        samtools \
            index \
            -b \
            -@ ${task.cpu} \
            ${sorted_sam}
        """

}


// process collapse {

//     tag "cDNA_Cupcake Collapse"
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.gff", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.gff" : "hq_unpolished_mapped_sorted_collapsed.gff" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.gff.unfuzzy", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.gff.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.gff.unfuzzy" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.group.txt" : "hq_unpolished_mapped_sorted_collapsed.group.txt" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt.unfuzzy", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.group.txt.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.group.txt.unfuzzy" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.rep.fa", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.rep.fa" : "hq_unpolished_mapped_sorted_collapsed.rep.fa" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.ignored_ids.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.ignored_ids.txt" : "hq_unpolished_mapped_sorted_collapsed.ignored_ids.txt" }

//     input:
//         tuple val(sample_name), file(sorted_sam), file(hq_reads)                    from sorted_for_collapsing_bam_ch
//         // file filtered from filtered_ch
//         // file fasta from unpolished_hq_fa_for_cupcake_ch

//     output:
//         tuple val(sample_name), file("*.collapsed.rep.fa"), file("*.collapsed.gff") into collapsed_seq_ch
//         file("*.unfuzzy")                                                           into collapsed_unfuzzy_ch
//         tuple val(sample_name), file("*.group.txt")                                 into collapsed_group_ch
//         file("*.group.txt.unfuzzy")                                                 into collapsed_group_unfuzzy_ch
//         file("*.ignored_ids.txt")                                                   into collapsed_ignored_ids_ch
    
//     script:
//         """
//         collapse_isoforms_by_sam \
//             --input ${hq_reads} \
//             --sam ${sorted_sam} \
//             --dun-merge-5-shorter \
//             --prefix ${hq_reads.baseName}
//         """
// }


// process abundance {

//     tag "Post-collapse abundance"

//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*abundance.txt", overwrite: true, { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.abundance.txt" : "hq_unpolished_mapped_sorted_collapsed.abundance.txt" }
//     publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*read_stat.txt", overwrite: true, { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed.read_stat.txt" : "hq_unpolished_mapped_sorted_collapsed.read_stat.txt" }

//     input:
//         tuple val(sample_name), file(group_file), file(cluster_report) from collapsed_group_ch.join(clustered_report_ch)

//     output:
//         tuple val(sample_name), file("*.abundance.txt")                into abundance_ch
//         file("*.read_stat.txt")                                        into read_stat_ch

//     script:
//     """
//         get_abundance_post_collapse ${group_file} ${cluster_report}
//     """

// }


// process filter_degraded {
//     // run filter_away_subset to remove 5' degraded transcripts

//     tag "Filter bad 5p ends"

//     input:
//         tuple val(sample_name), file(rep_file), file(gff_file), file(count_file) from collapsed_seq_ch.join(abundance_ch)

//     output:
//         file("*.filtered.abundance.txt")                                         into filtered_abundance_ch
//         file("*.filtered.gff")                                                   into filtered_annotation_ch
//         tuple val(sample_name), file("*.filtered.rep.fa")                        into filtered_seq_ch

//     script:
//         """
//         filter_away_subset ${count_file} ${gff_file} ${rep_file}
//         """

// }


// process sqanti {

//     tag "SQANTI3"
//     // errorStrategy 'ignore' // sometimes, there's just this one file...

//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.params.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep.params.txt": "hq_unpolished_mapped_sorted_collapsed_rep.params.txt" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep.renamed.fasta", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_renamed.fasta": "hq_unpolished_mapped_sorted_collapsed_rep_renamed.fasta" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_classification.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_classification.txt": "hq_unpolished_mapped_sorted_collapsed_rep_classification.txt" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.faa", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.faa": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.faa" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.fnn", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fnn": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fnn" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.fasta", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fasta": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fasta" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.fasta.fai", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.fasta.fai": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.fasta.fai" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.genePred", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.genePred": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.genePred" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.gtf", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.gtf": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.gtf" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.cds.gff", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.cds.gff": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.cds.gff" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected.sam", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected.sam": "hq_unpolished_mapped_sorted_collapsed_rep_corrected.sam" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_corrected_indels.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_corrected_indels.txt": "hq_unpolished_mapped_sorted_collapsed_rep_corrected_indels.txt" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep_junctions.txt", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep_junctions.txt": "hq_unpolished_mapped_sorted_collapsed_rep_junctions.txt" }
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "*.rep.genePred", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep.genePred": "hq_unpolished_mapped_sorted_collapsed_rep.genePred" }
//     publishDir path: "${params.sqanti}/${sample_name}/report_plot_data", mode: "copy", pattern: "*.csv", overwrite: true
//     publishDir path: "${params.sqanti}/${sample_name}",                  mode: "copy", pattern: "report.pdf", overwrite: true

//     input:
//         tuple val(sample_name), file(fixed_name_fa) from filtered_seq_ch

//     // Until we have a next step, these just keep the parent process from completing
//     output:
//         file("*.params.txt")                                     into sqanti_ch_1
//         file("*.rep.renamed.fasta")                              into sqanti_ch_2
//         tuple val(sample_name), file("*.rep_classification.txt") into sqanti_classification_ch
//         file("*.rep_corrected.faa")                              into sqanti_ch_4
//         file("*.rep_corrected.fnn")                              into sqanti_ch_5
//         file("*.rep_corrected.fasta")                            into sqanti_ch_6
//         file("*.rep_corrected.fasta.fai")                        into sqanti_ch_7
//         file("*.rep_corrected.genePred")                         into sqanti_ch_8
//         file("*.rep_corrected.gtf")                              into sqanti_corrected_annotation_ch
//         file("*.rep_corrected.cds.gff")                          into sqanti_ch_10
//         file("*.rep_corrected.sam")                              into sqanti_ch_11
//         file("*.rep_corrected_indels.txt")                       into sqanti_ch_12
//         file("*.rep_junctions.txt")                              into sqanti_junctions_ch
//         file("*.rep.genePred")                                   into sqanti_ch_14
//         file("*.csv")                                            into sqanti_plot_data_ch
//         file("report.pdf")                                       into sq_report_ch

//     script:
//       """
//       sqanti3_qc \
//           ${fixed_name_fa} \
//           ${params.annotation} \
//           ${params.genome} \
//           --cage_peak ${params.cage_peaks} \
//           --polyA_motif_list ${params.polyA_list} \
//           --cpus ${task.cpus} \
//           -vvv
//       """

// }


// process isoannot {

//     tag "annotation"

//     publishDir path: "${params.isoannot}/${sample_name}", mode: "copy", pattern: "*.gtf", overwrite: true, saveAs: { filename -> params.polish ? "hq_polished_mapped_sorted_collapsed_rep.isoannot.gtf": "hq_unpolished_mapped_sorted_collapsed_rep.isoannot.gtf" }

//     input:
//         file(corrected_annotation) from sqanti_corrected_annotation_ch
//         tuple val(sample_name), file(classification) from sqanti_classification_ch
//         file(junctions) from sqanti_junctions_ch

//     output:
//         file("*.gtf") into isoannot_annotation_ch

//     script:
//         """
//         IsoAnnotLite ${corrected_annotation} ${classification} ${junctions} --gff3 ${params.tappas_annotation} --output ${sample_name}.gtf
//         """

// }


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
