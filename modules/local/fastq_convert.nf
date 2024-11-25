process FASTQ_CONVERT{
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::pysam==0.22.1'
    container 'biocontainers/pysam==0.22.1'

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}/${meta.id}_bulk_m6A_R1.fastq.gz"), emit:convert_R1
    tuple val(meta), path("${meta.id}/${meta.id}_bulk_m6A_RC_R2.fastq.gz"), emit:convert_R2

    script:
    def prefix = "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    """
    fastq_convert.py \\
        --sample $prefix \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
    """
}