def run_draw_dataset(snakemake):
    """
    given an input probabilistic dataset, draw hard call
    output according to the input probability distributions
    """
    input_filename = snakemake.input
    input_format = snakemake.params["input_format"]
    output_filename = snakemake.output
    output_format = snakemake.params["output_format"]

    # for the moment, just touch output
    with open(output_filename, "w") as f:
        if input_filename or input_format or output_format or f:
            pass


try:
    snakemake
except NameError:
    pass
else:
    run_draw_dataset(snakemake)  # noqa: F821
