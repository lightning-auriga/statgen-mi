rule draw_from_minimac4:
    """
    create drawn dataset files from minimac4 data
    """
    input:
        lambda wildcards: config["imputed_datasets"][wildcards.dataset]["filename"][
            int(wildcards.filecode)
        ]
        if isinstance(config["imputed_datasets"][wildcards.dataset]["filename"], list)
        else config["imputed_datasets"][wildcards.dataset]["filename"],
    output:
        temp(
            "results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/data.vcf"
        ),
    params:
        input_format="minimac4",
        output_format="vcf",
    threads: 1
    resources:
        time="1:00:00",
        mem_mb="8000M",
        partition=config["queue"]["large_partition"],
    script:
        "../scripts/draw_dataset.py"
