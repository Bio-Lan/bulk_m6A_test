process BAM_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools==1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam_sorted)

    output:
    tuple val(meta), path("${meta.id}/*.bam"), emit:well_bam
    path "versions.yml", emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}/%!.bam"
    """
    mkdir ${meta.id}
    samtools \\
        split \\
        -d CB \\
        -M 400 \\
        -f $prefix \\
        --threads $task.cpus \\
        $bam_sorted
    
    rm ${meta.id}/-.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}