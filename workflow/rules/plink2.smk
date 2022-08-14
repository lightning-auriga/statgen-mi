rule multiple_imputation_plink2:
    """
    combine the results of plink2 analysis based on
    multiple imputation theory
    """
    input:
        lambda wildcards: expand(
            "results/{{analysis}}/{{dataset}}/{{tool}}/mi_runs/{runnum}/results.{{pheno_name}}.glm.{suffix}",
            runnum=[i + 1 for i in range(config["tools"][wildcards.tool]["mi_draws"])],
        ),
    output:
        "results/{analysis}/{dataset}/{tool}/output.tsv",
    params:
        tool=lambda wildcards: wildcards.tool,
        model=lambda wildcards: wildcards.suffix,
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
        vcf="results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/data.vcf.gz",
    output:
        filename="results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results.{pheno_name}.glm.linear",
        log=temp(
            "results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results.{pheno_name}.log"
        ),
    conda:
        "../envs/plink2.yaml"
    params:
        plink2_memlimit=config["plink2"]["memlimit"],
        plink2_outprefix="results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results",
        plink2_cols="chrom,pos,ref,alt,a1freq,a1count,test,nobs,orbeta,se,p",
        plink2_regression_modifiers="",
    threads: config["plink2"]["maxthreads"]
    resources:
        time="1:00:00",
        mem_mb=str(config["plink2"]["memlimit"]) + "M",
        partition=config["queue"]["large_partition"],
    shell:
        "plink2 --memory {params.plink2_memlimit} --threads {threads} "
        "--vcf {input.vcf} "
        "--glm hide-covar cols={params.plink2_cols} {params.plink2_regression_modifiers} "
        "--out {params.plink2_outprefix}"


use rule run_plink2_linear_regression as run_plink2_logistic_regression with:
    output:
        filename="results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results.{pheno_name}.glm.logistic.hybrid",
        log=temp(
            "results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results.{pheno_name}.log"
        ),
    params:
        memlimit=config["tools"]["plink"]["memlimit"] / 4,
        plink2_outprefix="results/{analysis}/{dataset}/{tool}/mi_runs/{runnum}/results",
        plink2_cols="chrom,pos,ref,alt,a1freq,a1freqcc,a1count,a1countcc,test,nobs,orbeta,se,p",
        plink2_regression_modifiers="--1",
