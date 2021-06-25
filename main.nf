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
      --ccs                 Default: project/02_ccs
      --demux               Default: project/03_demux
      --refined             Default: project/04_refined
      --clustered           Default: project/05_clustered
      --unpolished          Default: project/06_unpolished
      --polished            Default: project/07_polished
      --mapped              Default: project/08_mapped
      --sorted              Default: project/09_sorted
      --transcompressed     Default: project/10_transcompress
      --filtered            Default: project/11_filtered
      --collapsed           Default: project/12_collapsed
      --sqanti              Default: project/13_sqanti
      --isoannot            Default: project/14_isoannot

      NOTE: currently, chaining the results from cDNA_Cupcake has to be performed manually
    
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

params.input = "${params.raw}/*.subreads.bam"

Channel
    .fromPath( params.input, checkIfExists: true )
    .into{ raw_subreads_1; raw_subreads_2; raw_subreads_3 }

Channel
    .fromPath( params.barcodes, checkIfExists: true )
    .set{ barcodes_ch }


// CCS generation
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
        file("*.bam")            into ccs_calling_chunks_ch
        file("*.bam")            into ccs_calling_chunks_index_ch
        file("*.ccs_report.txt") into ccs_calling_report_ch
        file("*.json.gz")        into ccs_calling_metrics_ch
        file("*.log")            into ccs_calling_log_ch
        file(index)              into raw_index_ch_2

    script:
        if (params.runid != "")
            """
            ccs \
                --skip-polish \
                --min-passes=0 \
                --min-snr 4 \
                --min-rq 0.8 \
                --report-file ccs_report.${chunk}.txt \
                --report-json ccs_report.${chunk}.json \
                --metrics-json ccs_metrics.${chunk}.json \
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
        file(bam_chunks)  from ccs_calling_chunks_ch.collect()

    output:
        file("*.ccs.bam") into ccs_ch
        file("*.pbi")     into index_ch //delete?

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


// Classify
process demux { 

    tag "Demultiplexing samples"

    publishDir path: "${params.demux}", mode: "copy", pattern: "*.bam", overwrite: true
    publishDir path: "${params.logs}/demux", mode: "copy", pattern: "*.log", overwrite: true
    
    input:
        // val sample from demux_sample_name_ch
        file(called_ccs) from ccs_ch
        file(barcodes)   from barcodes_ch

    output:
        file("*.bam")                  into demuxed_bam_ch
        file("*.bam.pbi")              into demuxed_bam_index_ch
        file("*.consensusreadset.xml") into demuxed_readset_ch
        file("*.log")                  into demuxed_log_ch

    script:
        if (params.biosample != "")
            """
            lima \
                --isoseq \
                --peek \
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
                --peek \
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
        file(demuxed) from demuxed_bam_ch.flatten()
        // file(barcodes) from refine_barcodes_ch
    
    output:
        tuple val(sample_name), file("*.bam") into refined_ch
        file("*.pbi") into refined_index_ch
        file("*.csv") into refined_report_ch
        file("*.log") into refined_log_ch
        file("*.xml") into refined_consensusreadset_ch
        file("*.json") into refined_filter_ch

    script:
        sample_name = (demuxed.baseName =~ /demuxed\.([\s\S]*)_5p--[\s\S]*_3p/)[0][1]
        """
        isoseq3 \
            refine \
            --require-polya \
            --num-threads ${task.cpus} \
            --log-file refining.${sample_name}.log \
            --log-level ${params.refining_log_level} \
            --verbose \
            ${demuxed} \
            ${params.barcodes} \
            demuxed.refined.${sample_name}.flnc.bam
        """

}


// Cluster
process cluster {

    tag "Clustering"

    publishDir path: "${params.clustered}/${sample_name}", mode: "copy", pattern: "*.hq.bam", overwrite: true, saveAs: { filename -> "flnc.unpolished.hq.bam" }
    publishDir path: "${params.clustered}/${sample_name}", mode: "copy", pattern: "*.lq.bam", overwrite: true, saveAs: { filename -> "flnc.unpolished.lq.bam" }
    publishDir path: "${params.clustered}/${sample_name}", mode: "copy", pattern: "*.fasta.gz", overwrite: true, saveAs: { filename -> "flnc.unpolished.fasta.gz" }
    publishDir path: "${params.clustered}/${sample_name}", mode: "copy", pattern: "*.csv", overwrite: true, saveAs: { filename -> "flnc.unpolished.cluster_report.csv" }
    publishDir path: "${params.clustered}/${sample_name}", mode: "copy", pattern: "*.cluster", overwrite: true, saveAs: { filename -> "flnc.unpolished.cluster" }
    publishDir path: "${params.logs}/clustered", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "clustering.${sample_name}.log" }

    input:
        tuple val(sample_name), file(refined_bam)                                  from refined_ch
    
    output:
        tuple val(sample_name), file("*.hq.unpolished.bam"), file("*.hq.fasta.gz") into clustered_hq_results_ch
        tuple val(sample_name), file("*.hq.unpolished.bam"), file("*.hq.fasta.gz") into clustered_hq_results_skip_polish_ch
        file("*.unpolished.hq.bam.pbi")                                            into clustered_hq_index_ch
        tuple file("*.unpolished.lq.bam"), file("*.unpolished.lq.bam.pbi")         into clustered_lq_results_ch
        file("*.lq.fasta.gz")                                                      into clustered_lq_reads_ch
        file("*.cluster")                                                          into clustered_clusters_ch
        file("*.cluster_report.csv")                                               into clustered_report_ch
        file("*.log")                                                              into clustered_log_ch
        file("*.transcriptset.xml")                                                into clustered_transcriptset_ch

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
    publishDir path: "${params.logs}/polished", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "polish.${sample_name}.log" }

    input:
        tuple val(sample_name), file(unpolished_bam), file(hq_reads) from clustered_hq_results_ch
        file(index)                                                  from raw_index_ch_2
        file(subreads)                                               from raw_subreads_3
    
    output:
        tuple val(sample_name), file("*.bam"), file(hq_reads)        into polished_bam_ch
        file("*.log")                                                into polished_log_ch

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
                ${subreads} \
                ${unpolished_bam.baseName}.polished.bam
            """
        else
            """
            mv ${unpolished_bam} demuxed.unpolished.bam
            touch ${unpolished_bam.baseName}.log
            """

}


process mapping {

    tag "Mapping"

    publishDir path: "${params.mapped}/${sample_name}", mode: "copy", pattern: "*.bam", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped.sam" : "hq_unpolished_mapped.bam" }
    publishDir path: "${params.logs}/mapped", mode: "copy", pattern: "*.log", overwrite: true, saveAs: { filename -> "mapping.${sample_name}.log" }

    input:
        tuple val(sample_name), file(polished_bam), file(hq_reads) from polished_bam_ch
        tuple val(sample_name_alt), file(unpolished_bam)           from clustered_hq_results_skip_polish_ch

    output:
        tuple val(sample_name), file("*.bam"), file(hq_reads)      into mapped_bam_ch
        file("*.bai")                                              into mapped_index_ch
        file("*.alignment.log")                                    into mapped_log_ch

    script:
        if (params.polish)
            """
            pbmm2 \
                align \
                    --sample ${sample_name} \
                    --num-threads ${task.cpus} \
                    --preset ISOSEQ \
                    --log-level ${params.alignment_log_level} \
                    --log-file ${sample_name}.alignment.log \
                    --sort \
                    --bam-index BAI \
                    ${params.pbmm2_index} \
                    ${polished_bam} \
                    ${sample_name}.mapped.polished.flnc.bam
            """
        else
            """
            pbmm2 \
                align \
                    --sort \
                    --sort-threads ${task.cpus} \
                    --preset ISOSEQ \
                    --sample ${sample_name} \
                    --log-level ${params.alignment_log_level} \
                    --log-file ${sample_name_alt}.alignment.log \
                    --bam-index BAI \
                    ${params.pbmm2_index} \
                    ${unpolished_bam} \
                    ${sample_name_alt}.mapped.polished.flnc.bam
            """

    // script:
    //     """
    //     minimap2 \
    //         -H \
    //         -t ${task.cpus} \
    //         -ax splice:hq \
    //         -uf \
    //         --secondary=no \
    //         -O6,24 \
    //         -B4 \
    //         ${params.minimap_index} \
    //         ${polished_hq_fasta} \
    //         > ${polished_hq_fasta.baseName}.mapped.sam \
    //         2> ${polished_hq_fasta.baseName}.mapped.log
    //     """

}


process bam_to_sam {

    tag "BAM file conversion"

    input:
        tuple val(sample_name), file(mapped_bam), file(hq_reads) from mapped_bam_ch

    output:
        tuple val(sample_name), file("*.sam"), file(hq_reads)    into mapped_sam_ch

    script:
        if (params.ref_seq != "")
            """
            samtools \
                --threads \
                --write-index \
                --output-fmt SAM
                -o ${sample_name}.mapped.sam \
                --reference {params.ref_seq} \
                ${mapped_bam}
            """
        else
            """
            samtools \
                --threads \
                --write-index \
                --output-fmt SAM
                -o ${sample_name}.mapped.sam \
                ${mapped_bam}
            """

}


process collapse {

    tag "cDNA_Cupcake Collapse"
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.gff", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.gff" : "hq_unpolished_mapped_sorted_collapsed.gff" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.unfuzzy", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.unfuzzy" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.group.txt" : "hq_unpolished_mapped_sorted_collapsed.group.txt" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.group.txt.unfuzzy", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.group.txt.unfuzzy" : "hq_unpolished_mapped_sorted_collapsed.group.txt.unfuzzy" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.rep.fa", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.rep.fa" : "hq_unpolished_mapped_sorted_collapsed.rep.fa" }
    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.ignored_ids.txt", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed.ignored_ids.txt" : "hq_unpolished_mapped_sorted_collapsed.ignored_ids.txt" }

    input:
        tuple val(sample_name), file(sorted_sam), file(hq_reads) from mapped_sam_ch
        // file filtered from filtered_ch
        // file fasta from unpolished_hq_fa_for_cupcake_ch

    output:
        file("*.collapsed.gff")                                  into collapsed_annotation_ch
        tuple val(sample_name), file("*.collapsed.rep.fa")       into collapsed_seq_ch
        file("*.unfuzzy")                                        into collapsed_unfuzzy_ch
        tuple val(sample_name), file("*.group.txt")              into collapsed_group_ch
        file("*.group.txt.unfuzzy")                              into collapsed_group_unfuzzy_ch
        file("*.ignored_ids.txt")                                into collapsed_ignored_ids_ch
    
    script:
        """
        collapse_isoforms_by_sam \
            --input ${hq_reads} \
            --sam ${sorted_sam} \
            --dun-merge-5-shorter \
            --prefix ${hq_reads.baseName}
        """
}


process abundance {

    tag "Post-collapse abundance"

    publishDir path: "${params.collapsed}/${sample_name}", mode: "copy", pattern: "*.txt", overwrite: true
    
    input:
        file(cluster_report)                     from clustered_report_ch
        tuple val(sample_name), file(group_file) from collapsed_group_ch

    output:
        file("*.abundance.txt")                  into abundance_ch
        file("*.read_stat.txt")                  into read_stat_ch

    script:
    """
        get_abundance_post_collapse ${group_file} ${cluster_report}
    """

}


process filter_degraded {
    // run filter_away_subset to remove 5' degraded transcripts

    tag "Filter bad 5' ends"

    input:
        file(count_file) from abundance_ch
        file(gff_file) from collapsed_annotation_ch
        tuple val(sample_name), file(rep_file) from collapsed_seq_ch

    output:
        file("*.filtered.abundance.txt") into filtered_abundance_ch
        file("*.filtered.gff") into filtered_annotation_ch
        tuple val(sample_name), file("*.filtered.rep.fa") into filtered_seq_ch

    script:
        """
        filter_away_subset ${count_file} ${gff_file} ${rep_file}
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
        tuple val(sample_name), file(fixed_name_fa) from filtered_seq_ch

    // Until we have a next step, these just keep the parent process from completing
    output:
        file("*.params.txt")               into sqanti_ch_1
        file("*.rep.renamed.fasta")        into sqanti_ch_2
        tuple val(sample_name), file("*.rep_classification.txt") into sqanti_classification_ch
        file("*.rep_corrected.faa")        into sqanti_ch_4
        file("*.rep_corrected.fnn")        into sqanti_ch_5
        file("*.rep_corrected.fasta")      into sqanti_ch_6
        file("*.rep_corrected.fasta.fai")  into sqanti_ch_7
        file("*.rep_corrected.genePred")   into sqanti_ch_8
        file("*.rep_corrected.gtf")        into sqanti_corrected_annotation_ch
        file("*.rep_corrected.cds.gff")    into sqanti_ch_10
        file("*.rep_corrected.sam")        into sqanti_ch_11
        file("*.rep_corrected_indels.txt") into sqanti_ch_12
        file("*.rep_junctions.txt")        into sqanti_junctions_ch
        file("*.rep.genePred")             into sqanti_ch_14
        file("report.pdf") into sq_report_ch

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


process isoannot {

    tag "annotation"

    publishDir path: "${params.isoannot}/${sample_name}", mode: "copy", pattern: "*.gtf", overwrite: true, saveAs: { filename -> params.polished ? "hq_polished_mapped_sorted_collapsed_rep.isoannot.gtf": "hq_unpolished_mapped_sorted_collapsed_rep.isoannot.gtf" }

    input:
        file(corrected_annotation) from sqanti_corrected_annotation_ch
        tuple val(sample_name), file(classification) from sqanti_classification_ch
        file(junctions) from sqanti_junctions_ch

    output:
        file("*.gtf") into isoannot_annotation_ch

    script:
        """
        IsoAnnotLite ${corrected_annotation} ${classification} ${junctions} --gff3 ${params.tappas_annotation} --output ${sample_name}.gtf
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
