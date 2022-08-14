def run_construct_trait_file(snakemake):
    """
    build a plink-style phenotype output file
    given input

    this doesn't probably do very much at the moment
    """
    input_filename = (snakemake.input[0],)
    output_filename = snakemake.output[0]
    phenotype = snakemake.params["phenotype"]
    covariates = snakemake.params["covariates"]

    # for now, just touch output
    with open(output_filename, "w") as f:
        if input_filename or phenotype or covariates or f:
            pass


try:
    snakemake
except NameError:
    pass
else:
    run_construct_trait_file(snakemake)  # noqa: F821
