// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CSS_CALLING {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::pbccs=6.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbccs:6.0.0--h9ee0642_2"
    } else {
        container "quay.io/biocontainers/pbccs:6.0.0--h9ee0642_2"
    }

    input:
        tuple val(meta), path(subreads)

    output:
        tuple val(meta), path("*.ccs.${chunk}.bam"), emit: chunks
        tuple val(meta), path("ccs.log"), emit: log

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
        """
        ccs \
            --min-rq 0.9 \
            --log-level DEBUG \
            --log-file ccs.log \
            --num-threads ${task.cpus} \
            --chunk ${chunk}/${chunks} \
            ${subreads} \
            ${params.runid}.ccs.${chunk}.bam
        ccs --version | head -1 | sed -z 's/[^"]*\([[:digit:]]\.[[:digit:]]\.[[:digit:]]\).*/\1/g' > ${software}.version.txt
        """
}