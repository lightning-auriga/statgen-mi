rule combine_multiple_imputation:
    """
    for imputed datasets with multiple input files,
    merge results files into a single file; or just
    copy the file over really
    """
    input:
        lambda wildcards: expand(
            "results/{{analysis}}/{{dataset}}/{{tool}}/{{model}}/{filecode}/{{pheno_name}}.{{suffix}}.tsv",
            filecode=[
                i
                for i in range(
                    len(config["imputed_datasets"][wildcards.dataset]["filename"])
                )
            ]
            if isinstance(
                config["imputed_datasets"][wildcards.dataset]["filename"], list
            )
            else 0,
        ),
    output:
        "results/{analysis}/{dataset}/{tool}/{model}/{pheno_name}.{suffix}.tsv",
    threads: 1
    resources:
        time="1:00:00",
        mem_mb="1000M",
        partition=config["queue"]["small_partition"],
    shell:
        "cat {input} | awk 'NR == 1 || ! /effect_allele/' > {output}"


rule multiple_imputation_plink2:
    """
    combine the results of plink2 analysis based on
    multiple imputation theory
    """
    input:
        lambda wildcards: expand(
            "results/{{analysis}}/{{dataset}}/{{tool}}/{{model}}/{{filecode}}/mi_runs/{runnum}/results.{{pheno_name}}.{{suffix}}",
            runnum=[i + 1 for i in range(config["tools"][wildcards.tool]["mi_draws"])],
        ),
    output:
        "results/{analysis}/{dataset}/{tool}/{model}/{filecode}/{pheno_name}.{suffix}.tsv",
    params:
        tool=lambda wildcards: wildcards.tool,
        model=lambda wildcards: wildcards.suffix,
    wildcard_constraints:
        suffix="glm.linear|glm.logistic.hybrid",
    threads: 1
    resources:
        time="1:00:00",
        mem_mb="8000M",
        partition=config["queue"]["large_partition"],
    script:
        "../scripts/combine_analysis_output.py"


rule run_plink2_linear_regression:
    """
    deploy plink2 linear regression for drawn data
    """
    input:
        vcf="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/data.vcf.gz",
        tbi="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/data.vcf.gz.tbi",
        pheno="results/{analysis}/{dataset}/{tool}/{model}/{pheno_name}.pheno",
    output:
        filename="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results.{pheno_name}.glm.linear",
        log=temp(
            "results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results.{pheno_name}.log"
        ),
    conda:
        "../envs/plink2.yaml"
    params:
        plink2_exec=config["tools"]["plink2"]["executable"],
        plink2_memlimit=config["tools"]["plink2"]["maxmem"],
        plink2_outprefix="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results",
        plink2_cols="chrom,pos,ref,alt,a1freq,a1count,test,nobs,orbeta,se,p",
        plink2_regression_modifiers="",
        phenos=lambda wildcards: " --pheno-name {} ".format(wildcards.pheno_name),
        covars=lambda wildcards: " --covar-name {} ".format(
            " ".join(config["regression_models"][wildcards.model]["covariates"])
        )
        if "covariates" in config["regression_models"][wildcards.model]
        else "",
    threads: config["tools"]["plink2"]["maxthreads"]
    resources:
        time="1:00:00",
        mem_mb=str(config["tools"]["plink2"]["maxmem"]) + "M",
        partition=config["queue"]["large_partition"],
    shell:
        "{params.plink2_exec} --memory {params.plink2_memlimit} --threads {threads} "
        "--vcf {input.vcf} --double-id "
        "--glm hide-covar cols={params.plink2_cols} {params.plink2_regression_modifiers} "
        "--pheno {input.pheno} {params.phenos} {params.covars} "
        "--out {params.plink2_outprefix} && mv {params.plink2_outprefix}.log {output.log}"


use rule run_plink2_linear_regression as run_plink2_logistic_regression with:
    output:
        filename="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results.{pheno_name}.glm.logistic.hybrid",
        log=temp(
            "results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results.{pheno_name}.log"
        ),
    params:
        plink2_exec=config["tools"]["plink2"]["executable"],
        plink2_memlimit=config["tools"]["plink2"]["maxmem"] / 4,
        plink2_outprefix="results/{analysis}/{dataset}/{tool}/{model}/{filecode}/mi_runs/{runnum}/results",
        plink2_cols="chrom,pos,ref,alt,a1freq,a1freqcc,a1count,a1countcc,test,nobs,orbeta,se,p",
        plink2_regression_modifiers="--1",
        phenos=lambda wildcards: " --pheno-name {} ".format(wildcards.pheno_name),
        covars=lambda wildcards: " --covar-name {} ".format(
            " ".join(config["regression_models"][wildcards.model]["covariates"])
        )
        if "covariates" in config["regression_models"][wildcards.model]
        else "",
