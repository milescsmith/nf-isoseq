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

extract_bc = { item -> 
    item =~ /bc(\d+)\-\[FR]/
}

process ccs_indexing {
    conda "bioconda::pbbam"

    tag "CCS indexing"
    publishDir '${params.raw}', mode: 'copy', pattern : '*.bam.pbi', overwrite: false
    
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

    tag "CCS calling"
    publishDir '${params.ccs}', mode: 'copy', pattern: '*.ccs.${chunk}.bam', overwrite: true
    publishDir '${params.ccs}', mode: 'copy', pattern: "*.log", overwrite: true

    input:
        file(subreads) from raw_subreads_2
        file(index) from raw_index_ch
        each chunk from 1..params.chunks
        val chunks from params.chunks

    output:
        file "*.ccs.${chunk}.bam" into ccs_chunks_ch
        file "ccs.log" into css_log

    script:
        """
        ccs \
            --min-rq 0.9 \
            --log-level DEBUG \
            --log-file ccs.log \
            --num-threads ${task.cpus} \
            --chunk ${chunk}/${chunks} \
            ${subreads} \
            ${params.runid}.ccs.bam
        """

}

process ccs_chunk_merging {
    conda "bioconda::pbccs"

    tag "CCS chunk merging"
    publishDir '${params.ccs}', mode: 'copy', pattern: "*.ccs.bam", overwrite: true
    publishDir '${params.ccs}', mode: 'copy', pattern: "*.pbi", overwrite: true

    input:
        file bam_chunks from ccs_chunks_ch.collect()

    output:
        file "*.ccs.bam" into ccs_ch
        file "*.pbi" into index_ch

    script:
        """
        pbmerge -o ${params.runid}.ccs.bam ${bam_chunks}
        pbindex ${params.runid}.ccs.bam
        """
}

process demux {
    conda "bioconda::lima"

    tag "Demultiplexing samples"
    publishDir '${params.demux}', mode: 'copy', pattern: '*.bam', overwrite: true
    publishDir '${params.demux}', mode: 'copy', pattern: '*.log', overwrite: true

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
            --split-bam \
            ${called_ccs} \
            ${params.barcodes} \
            demuxed.bam
        """
}

process refine {
    conda "bioconda::isoseq3"

    tag "Refining"
    publishDir '${params.refine}', mode: 'copy', pattern: '*.flnc.bam', overwrite: true
    publishDir '${params.refine}', mode: 'copy', pattern: '*.log', overwrite: true

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
    conda "bioconda::isoseq3"

    tag "Clustering"
    publishDir '${params.unpolished}', mode: 'copy', pattern: '*.bam', overwrite: true
    publishDir '${params.unpolished}', mode: 'copy', pattern: '*.log', overwrite: true

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
    conda "bioconda::isoseq3"

    tag "Polishing"
    publishDir '${params.polished}', mode: 'copy', pattern: '*.bam', overwrite: true
    publishDir '${params.polished}', mode: 'copy', pattern: '*.log', overwrite: true

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
    conda "bioconda::gmap==2020.04.08"

    tag "Mapping"
    publishDir '${params.mapped}', mode: 'copy', pattern: '*.sam', overwrite: true
    publishDir '${params.mapped}', mode: 'copy', pattern: '*.log', overwrite: true

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
    conda "bioconda::samtools==1.10"

    tag "Sorting"
    publishDir '${params.sorted}', mode: 'copy', pattern: '*.bam', overwrite: true

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
    conda "bioconda::samtools==1.10"

    tag "Transcompressing"
    publishDir '${params.transcompressed}', mode: 'copy', pattern: '*.fasta.gz', overwrite: true

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
    tag "Filtering"
    publishDir '${params.filtered}', mode: 'copy', pattern: '*.bam', overwrite: true
    
    container "registry.gitlab.com/milothepsychic/filter_sam"

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
    publishDir '${params.collapsed}', mode: 'copy', pattern: '*.gff', overwrite: true

    container "milescsmith/cdna_cupcake"

    input:
        // val sample from sample_name_ch
        file filtered from filtered_ch
        file fasta from bgzipped_hq_ch_2

    output:
        file "${sample}.gff" into collapsed_ch
    
    script:
    """
    collapse_isoforms_by_sam.py \
        --input ${fasta} \
        --sam ${filtered} \
        --dun-merge-5-shorter \
        --prefix ${sample}
    """
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
