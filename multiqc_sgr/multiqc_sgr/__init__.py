def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "bulk_m6A/stats": {
            "fn": "*bulk_m6A.*.stats.json"
        },
        "bulk_m6A/motif": {
            "fn": "*bulk_m6A.conversion.motif.json"
        },
        "bulk_m6A/boxplot": {
            "fn": "*bulk_m6A.substitution.boxplot.json"
        },
        "bulk_m6A/barplot": {
            "fn": "*bulk_m6A.substitution.barplot.json"
        },
        "bulk_m6A/well_inf": {
            "fn": "*bulk_m6A.quant.well_inf.json"
        },
        "bulk_m6A/label_rate": {
            "fn": "*bulk_m6A.quant.labeled_rate.json"
        }
    }
    config.update_dict(config.sp, sgr_search_patterns)
