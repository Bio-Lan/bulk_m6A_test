process CONVERSION{
    tag "$meta.id"
    label 'process_high'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 bioconda::samtools==1.20'
    container "qaqlans/sgr-accura-2"

    input:
    tuple val(meta), path(well_bam)

    output:
    tuple val(meta), path("${meta.id}/*.PosTag.bam"), emit:conv_wellbam
    tuple val(meta), path("${meta.id}/*.PosTag.bam.bai"), emit:conv_bambai
    tuple val(meta), path("${meta.id}/*.conversion_loci.csv"), emit:conv_wellloci

    script:
    def prefix = "${meta.id}"
    """
    conversion.py \\
        --sample $prefix \\
        --wellBAM $well_bam \\
        --GTF ${params.gtf} \\
        --fastafile ${params.fasta} \\
        --conversion_type ${params.conversion_type} \\
        --basequalilty ${params.basequalilty} \\
        --MAPQ ${params.MAPQ}
    """
}