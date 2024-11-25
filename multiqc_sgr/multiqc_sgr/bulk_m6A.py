import json
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table, bargraph, box

ASSAY = "bulk_m6A"

# Initialise the logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name = ASSAY,
            anchor = ASSAY,
            info = "mapping, demultiplexing and quantification for bulk_dynaseq sequencing data",
        )
        log.info(f"Running module: {self.name}")

        stat_data = self.parse_json(self.name, "stats")
        motif_data = self.parse_json(self.name, "motif")
        box_data = self.parse_json(self.name, "boxplot")
        bar_data = self.parse_json(self.name, "barplot")
        labeled_data = self.parse_json(self.name, "label_rate")
        table_data = self.parse_json(self.name, "well_inf")
        
        if all(len(x) == 0 for x in [stat_data, box_data, bar_data, labeled_data, table_data]):
            raise ModuleNoSamplesFound
        
        sample_list = list(stat_data.keys())
        sample_list.sort()

        self.general_stats_table(stat_data)

        lable_helptext = """
        Only output wells that meet the conditions(default: UMI>=500,gene>=0 ).  
        If no well pass the filter,output wells that UMI and gene >= 0.
        """      

        table_helptext = """
        The rowname of the table is the well barcode sequence. Only output wells that meet the conditions(default: UMI>=500,gene>=0 ).  
        If no well pass the filter,output wells that UMI and gene >= 0.
        """      
        
        for sample in sample_list:
            self.add_section(
                name=f"{sample} - DRACH motifs",
                anchor=f"{sample}_motif",
                plot=self.motif_bar(motif_data[sample],sample)
            )

            self.add_section(
                name=f"{sample} - substitution rate per well",
                anchor=f"{sample}_barplot",
                plot=self.substitution_bar(bar_data[sample],sample)
            )
            self.add_section(
                name=f"{sample} - substitution rate per type",
                anchor=f"{sample}_boxplot",
                plot=self.substitution_box(box_data[sample],sample)
            )
            self.add_section(
                name=f"{sample} - labeled per well",
                anchor=f"{sample}_well_labeled",
                helptext=lable_helptext,
                plot=self.quant_labeled(labeled_data[sample],sample)
            )
            self.add_section(
                name=f"{sample} - per well",
                anchor=f"{sample}_well_table",
                helptext=table_helptext,
                plot=self.quant_table(table_data[sample],sample)
            )

    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict
    
    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green"
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Q30 Bases in CB+UMI": {
                "title": "Q30 Bases in CB+UMI",
                "description": "Q30 Bases in CB+UMI",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Q30 Bases in RNA read": {
                "title": "Q30 Bases in RNA read",
                "description": "Q30 Bases in RNA read",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Reads Mapped To Unique Loci": {
                "title": "Unique Reads",
                "description": "Percent of valid reads mapped to unique loci on the genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Reads Mapped To Multiple Loci": {
                "title": "Multi Reads",
                "description": "Percent of valid reads mapped to multiple loci on the genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Reads Mapped Uniquely To Transcriptome": {
                "title": "Counted Unique Reads",
                "description": "Percent of valid reads mapped uniquely to transcriptome; These reads are used for UMI counting",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green"
            },
            "Mapped Reads Assigned To Exonic Regions": {
                "title": "Exonic Reads",
                "description": "Percent of mapped reads assigned to exonic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned To Intronic Regions": {
                "title": "Intronic Reads",
                "description": "Percent of mapped reads assigned to intronic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned To Intergenic Regions": {
                "title": "Intergenic Reads",
                "description": "Percent of mapped reads assigned to intergenic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Mapped Reads Assigned Antisense To Gene": {
                "title": "Antisense Reads",
                "description": "Percent of mapped reads assigned antisense to gene",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True
            },
            "Median UMI across Well": {
                "title": "Median UMI across Well",
                "description": "The median number of UMI counts across wells.",
                "scale": "blue",
                "format": "{:,.0f}"
            },
            "Mean UMI across Well": {
                "title": "Mean UMI across Well",
                "description": "The mean number of UMI counts across wells.",
                "scale": "blue",
                "format": "{:,.0f}"
            },
            "Median Genes across Well": {
                "title": "Median Genes across Well",
                "description": "The median number of Genes counts across wells.",
                "scale": "blue",
                "format": "{:,.0f}"
            },
            "Mean Genes across Well": {
                "title": "Mean Genes across Well",
                "description": "The mean number of Genes counts across wells.",
                "scale": "blue",
                "format": "{:,.0f}"
            },
            "Median labeled rate across wells": {
                "title": "Median labeled rate across wells",
                "description": "The median UMI labeled rate across wells.",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green"
            },
            "Mean labeled rate across wells": {
                "title": "Mean labeled rate across wells",
                "description": "The mean UMI labeled rate across wells.",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green"
            }
        }
        self.general_stats_addcols(summary_data, headers=headers)
    
    def motif_bar(self,motif_data,sample):
        pconfig = {
            "id": f"{sample}-motif",
            "title": "DRACH motifs",
            "ylab": "Number of m6A site",
            "cpswitch": False,
            "sort_samples": False
        }
        return bargraph.plot(motif_data, pconfig=pconfig)
    
    def substitution_bar(self,barplot_data,sample):
        pconfig = {
            "id": f"{sample}-barplot",
            "title": "Barplot for Substitution per well",
            "ylab": "Nucleotide substitution rate",
            "cpswitch": False,
            "sort_samples":False
        }
        cats = [
            {
                "C_to_T": {'color': '#ae017e'},
                "A_to_C": {'color': '#c2e699'},
                "A_to_G": {'color': '#78c679'},
                "A_to_T": {'color': '#238443'},
                "C_to_A": {'color': '#fbb4b9'},
                "C_to_G": {'color': '#f768a1'},
                "G_to_A": {'color': '#bdd7e7'},
                "G_to_C": {'color': '#6baed6'},
                "G_to_T": {'color': '#2171b5'},
                "T_to_C": {'color': '#fb6a4a'},
                "T_to_A": {'color': '#fcae91'},
                "T_to_G": {'color': '#cb181d'}
            }
        ]
        return bargraph.plot(barplot_data, cats=cats, pconfig=pconfig)
    
    def substitution_box(self,boxplot_data,sample):
        pconfig = {
            "id": f"{sample}-boxplot",
            "title": "Boxplot for Substitution per well",
            "xlab": "Nucleotide substitution rate"
        }
        return box.plot(boxplot_data, pconfig=pconfig)

    def quant_labeled(self,label_data,sample):
        pconfig = {
            "id": f"{sample}-labeled",
            "title": "Labeled fraction per well",
            "xlab": "Fraction of labeled"
        }
        return box.plot(label_data,pconfig=pconfig)
    
    def quant_table(self,table_data,sample):
        table_config = {
            "id": f"{sample}-table",
            "title":"Detailed information per Well",
            "namespace": f"{sample}",
            "col1_header":"Well",
            "sort_rows":False
        }
        headers = {
            "UMI": {
                "title": "UMI",
                "format": "{:,.0f}"
            },
            "gene": {
                "title": "gene",
                "format": "{:,.0f}"
            },
            "labeled_UMI": {
                "title": "labeled_UMI",
                "format": "{:,.0f}"
            },
            "labeled_gene": {
                "title": "labeled_gene",
                "format": "{:,.0f}"
            }
        }
        return table.plot(table_data,pconfig=table_config,headers=headers)
    