import pandas as pd


def run_construct_trait_file(snakemake):
    """
    build a plink-style phenotype output file
    given input

    this doesn't probably do very much at the moment
    """
    input_filename = snakemake.input[0]
    output_filename = snakemake.output[0]
    phenotype = snakemake.params["phenotype"]
    covariates = snakemake.params["covariates"]
    data = pd.read_table(input_filename, sep="\t")
    if phenotype not in data.columns:
        raise ValueError("requested phenotype not present in model file")
    if covariates is not None:
        if not all(item in data.columns for item in covariates):
            raise ValueError("requested covariates not present in model file")
    data.to_csv(output_filename, sep="\t", index=False)


try:
    snakemake
except NameError:
    pass
else:
    run_construct_trait_file(snakemake)  # noqa: F821
