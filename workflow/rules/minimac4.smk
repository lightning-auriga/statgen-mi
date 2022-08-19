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
        "results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/data.vcf.gz",
    params:
        input_format="minimac4",
        output_format="vcf",
        bcftools_exec=config["tools"]["bcftools"]["executable"],
        bcftools_plugin_path=config["tools"]["bcftools"]["plugin_path"],
    threads: 1
    resources:
        time="1:00:00",
        mem_mb="8000M",
        partition=config["queue"]["large_partition"],
    shell:
        'BCFTOOLS_PLUGINS="{params.bcftools_plugin_path}" {params.bcftools_exec} +setGT {input} -O z -o {output} -- -t a -n r -r GP'
