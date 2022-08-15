rule bgzip_file:
    input:
        "{filename}",
    output:
        "{filename}.gz",
    conda:
        "../envs/tabix.yaml"
    threads: 1
    resources:
        time="0:30:00",
        mem_mb="1000M",
        partition=config["queue"]["small_partition"],
    shell:
        "bgzip -c {input} > {output}"


rule tabix_vcf:
    """
    standard rule for indexing compressed vcf
    """
    input:
        "{filename}.vcf.gz",
    output:
        "{filename}.vcf.gz.tbi",
    conda:
        "../envs/tabix.yaml"
    threads: 1
    resources:
        time="0:30:00",
        mem_mb="1000M",
        partition=config["queue"]["small_partition"],
    shell:
        "tabix -p vcf {input}"
