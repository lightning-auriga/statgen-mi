import math
import re
from contextlib import ExitStack


def plink2_handler(lines, model_name):
    """
    handle MI combination for results from plink2 glm
    """
    varid = ""
    effect_allele = ""
    sample_size = 0
    betas = []
    variances = []
    pattern = None
    if model_name == "linear":
        pattern = re.compile(
            "[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t"
            "[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t"
        )
    elif model_name == "logistic":
        pattern = re.compile(
            "[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t"
            "[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t"
            "([^\t]+)\t([^\t]+)\t"
        )
    else:
        raise ValueError(
            "combine_analysis_output plink2 handler does not currently support model {}".format(
                model_name
            )
        )
    for line in lines:
        match = pattern.match(line)
        if match:
            varid = match.group(1)
            effect_allele = match.group(2)
            if match.group(5) != "NA":
                sample_size = int(match.group(3))
                if model_name == "logistic":
                    betas.append(math.log(float(match.group(4))))
                else:
                    betas.append(float(match.group(4)))
                variances.append(pow(float(match.group(5)), 2))
        else:
            raise ValueError(
                "plink2 regression result does not match expected format: {}".format(
                    line
                )
            )
    if len(betas) > 0:
        overall_beta = sum(betas) / len(betas)
        within_variance = sum(variances) / len(variances)
        between_variance = sum(
            [pow(overall_beta - current_beta, 2) for current_beta in betas]
        ) / len(betas)
        total_variance = (
            within_variance + (len(betas) + 1) / len(betas) * between_variance
        )
        fraction_total_missing = (
            between_variance * (len(betas) + 1 / len(betas))
        ) / total_variance
        dof_old = (len(betas) - 1) / pow(fraction_total_missing, 2)
        dof_obs = (
            (sample_size - 1)
            / (sample_size + 1)
            * (sample_size - 2)
            * (1 - fraction_total_missing)
        )
        dof_adj = dof_old * dof_obs / (dof_old + dof_obs)
        tstat = overall_beta / math.sqrt(total_variance)
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            varid,
            effect_allele,
            overall_beta,
            within_variance,
            between_variance,
            total_variance,
            dof_adj,
            tstat,
        )
    else:
        return "{}\t{}\tNA\tNA\tNA\tNA\tNA\n".format(varid, effect_allele)


def run_combine_analysis_output(snakemake):
    """
    given the results of multiple passes through an analysis tool,
    combine the results with MI logic
    """
    input_filenames = snakemake.input
    output_filename = snakemake.output[0]
    tool_name = snakemake.params["tool"]
    model_name = snakemake.params["model"]

    with ExitStack() as stack, open(output_filename, "w") as out:
        files = [
            stack.enter_context(open(filename, "r")) for filename in input_filenames
        ]
        # chomp header
        lines = [f.readline() for f in files]
        while True:
            lines = [f.readline() for f in files]
            if all(entry == "" for entry in lines):
                break
            elif any(entry == "" for entry in lines):
                raise ValueError("input files for combine_analysis_output are jagged")
            if tool_name == "plink2":
                formatted_line = plink2_handler(lines, model_name)
            else:
                raise ValueError(
                    "combine_analysis_output does not currently support tool {}".format(
                        tool_name
                    )
                )
            out.write(formatted_line)


try:
    snakemake
except NameError:
    pass
else:
    run_combine_analysis_output(snakemake)  # noqa: F821
