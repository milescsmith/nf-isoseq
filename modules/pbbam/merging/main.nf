// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PBBAM_MERGING {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::pbbam=1.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbbam:1.6.0--h058f120_1"
    } else {
        container "quay.io/biocontainers/pbbam:1.6.0--h058f120_1"
    }

    input:
        tuple val(meta), path(chunks)

    output:
        tuple val(meta), path("*.bam.pbi"), emit: index
        file "*.ccs.bam" into ccs_ch
        file "*.pbi" into index_ch

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
        """
        pbmerge -o ${params.runid}.ccs.bam ${chunks}
        pbmerge --version | sed -e 's/pbmerge //g' > ${software}.version.txt
        """
}