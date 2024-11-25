process SUBSTITUTION{
    tag "$meta.id"
    label 'process_medium'
    
    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 bioconda::samtools==1.20'
    container "qaqlans/sgr-accura-2"
    
    input:
    tuple val(meta), path(conv_sample), path(conv_wellbam,stageAs: "bam_file/*")

    output:
    tuple val(meta), path("*.csv"), emit: substitution_stat

    script:
    def prefix = "${meta.id}"
    def conv_wellbam = conv_wellbam.join(",")
    """
    substitution.py \\
        --sample $prefix \\
        --conv_csv $conv_sample \\
        --conv_bam $conv_wellbam
    """
}