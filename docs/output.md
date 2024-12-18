- [Main Output](#main-output)
- [modules](#modules)
  - [filter\_gtf](#filter_gtf)
  - [star\_genome](#star_genome)
  - [protocol\_cmd](#protocol_cmd)
  - [starsolo](#starsolo)
  - [conversion](#conversion)
  - [conversion\_summary](#conversion_summary)
  - [substitution](#substitution)
  - [quant](#quant)
  - [report\_summary](#report_summary)
  - [multiqc-sgr](#multiqc-sgr)
  - [pipeline\_info](#pipeline_info)
  - [fastqc(Optional)](#fastqcoptional)

# Main Output

- `multiqc/multiqc_report.html` HTML report containing QC metrics.
- `quant/{sample}.matrix` The expression matrix of total/labeled/unlabeled in Matrix Market Exchange Formats.
- `starsolo/{sample}.Aligned.sortedByCoord.out.bam` Bam file contains coordinate-sorted reads aligned to the genome.

# modules

## filter_gtf

This module has the same functionality as [`cellranger mkgtf`](https://kb.10xgenomics.com/hc/en-us/articles/360002541171-What-criteria-should-I-use-with-the-mkgtf-tool-when-making-a-custom-reference-for-Cell-Ranger)

> GTF files can contain entries for non-polyA transcripts that overlap with protein-coding gene models. These entries can cause reads to be flagged as mapped to multiple genes (multi-mapped) because of the overlapping annotations. In the case where reads are flagged as multi-mapped, they are not counted.

> We recommend filtering the GTF file so that it contains only gene categories of interest by using the cellranger mkgtf tool. Which genes to filter depends on your research question.

The filtering criteria is controlled by the argument `--keep_attributes`. The default value of this argument is the same as the [reference used by cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_3.0.0)

> [!NOTE]
> gtf files from [genecode](https://www.gencodegenes.org/) use `gene_type` instead of `gene_biotype`.
>
> ```
> --keep_attributes "gene_type=protein_coding,lncRNA..."
> ```

**Output files**

- `*.filtered.gtf` GTF file after filtering.
- `gtf_filter.log` log file containing number of lines filtered in the original gtf file.



## star_genome

Generate STAR genome index. Detailed documents can be found in the [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

> [!TIP]
> Once you have the indices from a workflow run you should save them somewhere central and reuse them in subsequent runs using custom config files or command line parameters.

**Output files**

- `{genome_name}/` STAR genome index folder.



## protocol_cmd

Generate STARSolo command-line arguments accordingly.

**Output files**

- `{sample}.protocol_cmd.txt` STARSolo command-line arguments.
- `{sample}.bulk_m6A..protocol.stats.json` Protocol.



## starsolo

Descriptions of parameters and files can be found in [STARSolo documents](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

When you have questions, [STAR’s github issue](https://github.com/alexdobin/STAR/issues) is also a great place to find answers and help.

> [!NOTE]
> The command line arguments in this STARsolo documentation may not be up to date. For the latest STARSolo arguments, please refer to The STAR Manual.

**Output files**

- `{sample}.Aligned.sortedByCoord.out.bam` Bam file contains coordinate-sorted reads aligned to the genome.
- `{sample}.Solo.out` Starsolo output file.



## conversion

Get conversion pos in each read.

**Output files**

- `{sample}.PosTag.bam` Bam file with conversion info.
- `{well}.conversion_loci.csv` CT conversion sites info in csv format.



## conversion_summary

Filter conversion site

**Output files**
- `{well}.m6A_loci.csv` Sites were filtered and identified as m6A site.
- `{well}.loci_in_gene.csv` m6A site number in gene.
- `{well}.loci_motif.csv` m6A site DRACH motif.
- `{sample}.conversion.csv` The numbers of CT conversion sites and m6A sites of each well.



## substitution

Computes the overall conversion rates in reads.

**Output files**

- `{sample}.substitution.csv` Comma-separated table of the overall conversion rates.



## quant

Quantify total RNA, unlabeled and labeled RNA.

**Output files**

- `{sample}_{RNA}` The expression matrix of total/labeled/unlabeled in Matrix Market Exchange Formats.



## report_summary

Get multiqc report data.

**Output files**

- `{sample}.*.json` Summary file for pipeline analysis steps.


## multiqc-sgr

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

[multiqc-sgr](https://pypi.org/project/multiqc-sgr/) adds some modules on this basis to facilitate the visualization of single cell-related data.

**Output files**

- `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
- `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
- `multiqc_plots/`: directory containing static images from the report in various formats.



## pipeline_info

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files**

- Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
- Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
- Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
- Parameters used by the pipeline run: `params.json`.

## fastqc(Optional)

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files**

- `*_fastqc.html`: FastQC report containing quality metrics.
- `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
