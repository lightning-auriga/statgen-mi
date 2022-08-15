rule construct_phenotype_data:
    """
    build plink-format phenotype/covariate file for analysis
    """
    input:
        filename=lambda wildcards: config["regression_models"][wildcards.model][
            "filename"
        ],
    output:
        filename="results/{analysis}/{dataset}/{tool}/{model}/{pheno_name}.pheno",
    params:
        phenotype=lambda wildcards: config["regression_models"][wildcards.model][
            "phenotype"
        ],
        covariates=lambda wildcards: config["regression_models"][wildcards.model][
            "covariates"
        ]
        if "covariates" in config["regression_models"][wildcards.model]
        else None,
    threads: 1
    resources:
        time="1:00:00",
        mem_mb="1000M",
        partition=config["queue"]["small_partition"],
    script:
        "../scripts/construct_trait_file.py"
