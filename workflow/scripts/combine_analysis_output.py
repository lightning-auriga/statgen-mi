def run_combine_analysis_output(snakemake):
    """
    given the results of multiple passes through an analysis tool,
    combine the results with MI logic
    """
    input_filenames = snakemake.input
    output_filename = snakemake.output[0]
    tool_name = snakemake.params["tool"]
    model_name = snakemake.params["model"]

    # for the moment, just touch output
    with open(output_filename, "w") as f:
        for input_filename in input_filenames:
            if tool_name or model_name or f:
                pass


try:
    snakemake
except NameError:
    pass
else:
    run_combine_analysis_output(snakemake)  # noqa: F821
