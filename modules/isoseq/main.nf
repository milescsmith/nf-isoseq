// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LIMA_DEMUXING {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::lima==2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lima:2.0.0--0"
    } else {
        container "quay.io/biocontainers/lima:2.0.0--0"
    }

    input:
        tuple val(meta), path(called_css)
        tuple val(meta), val(barcodes)

    output:
        tuple val(meta), path("demuxed.bam"), emit: demuxed

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
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

        lima --version | sed -z 's/[^"]*\([[:digit:]]\.[[:digit:]]\.[[:digit:]]\).*/\1/g' > ${software}.version.txt
        """
}