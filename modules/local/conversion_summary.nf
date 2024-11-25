process CONVERSION_SUMMARY{
    tag "$meta.id"
    label 'process_median'

    conda 'conda-forge::pandas==2.2.1 bioconda::pysam==0.22.1 bioconda::samtools==1.20'
    container "qaqlans/sgr-accura-2"

    input:
    tuple val(meta), path(well_loci, stageAs: "?/*")
    
    output:
    tuple val(meta), path("${meta.id}/")
    tuple val(meta), path("${meta.id}/*.m6A_loci.csv"), emit: conv_loci
    tuple val(meta), path("${meta.id}/*.loci_motif.csv"), emit: conv_motif
    tuple val(meta), path("${meta.id}.conversion.csv"), emit: conv_sample
    tuple val(meta), path("${meta.id}.check_base.csv"), optional:true

    script:
    def well_loci = well_loci.join(",")
    def args = params.m6A_check ? "--check" : ''
    """       
    conversion_summary.py \\
        --sample ${meta.id} \\
        --well_loci $well_loci \\
        --min_m6A_depth ${params.min_m6A_depth} \\
        --min_coverage ${params.min_coverage} \\
        --min_m6A_threshold ${params.min_m6A_threshold} \\
        --max_m6A_threshold ${params.max_m6A_threshold} \\
        $args
    """
}