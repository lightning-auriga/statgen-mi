import math
import re
from contextlib import ExitStack


def format_mi_result(
    varid, effect_allele, other_allele, sample_size, betas, variances, mi_run_count
):
    """
    given parsed information from multiple association runs,
    combine and return MI results
    """
    output_float_dec = 6
    if len(betas) == mi_run_count:
        overall_beta = sum(betas) / len(betas)
        within_variance = sum(variances) / len(variances)
        between_variance = sum(
            [pow(overall_beta - current_beta, 2) for current_beta in betas]
        ) / len(betas)
        # set test statistic defaults assuming invariant draws, as will
        # happen for perfect confidence genotypes
        total_variance = within_variance
        dof_adj = 1
        tstat = overall_beta / math.sqrt(total_variance)
        if not math.isclose(between_variance, 0):
            total_variance = (
                within_variance + (len(betas) + 1) / len(betas) * between_variance
            )
            fraction_total_missing = (
                between_variance * (1 + 1 / len(betas))
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
        return "{}\n".format(
            "\t".join(
                [
                    varid,
                    effect_allele,
                    other_allele,
                    str(round(overall_beta, output_float_dec)),
                    str(round(within_variance, output_float_dec)),
                    str(round(between_variance, output_float_dec)),
                    str(round(total_variance, output_float_dec)),
                    str(round(dof_adj, output_float_dec)),
                    str(round(tstat, output_float_dec)),
                ]
            )
        )
    # for the time being, enforce all simulations converge. if they don't,
    # treat the result as not converged
    return "{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\n".format(varid, effect_allele, other_allele)


def plink2_handler(lines, model_name):
    """
    handle MI combination for results from plink2 glm
    """
    varid = ""
    ref_allele = ""
    alt_allele = ""
    effect_allele = None
    sample_size = 0
    betas = []
    variances = []
    pattern = None
    if model_name == "glm.linear":
        pattern = re.compile(
            "[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t"
            "[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t"
        )
    elif model_name == "glm.logistic.hybrid":
        pattern = re.compile(
            "[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t"
            "[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t"
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
            ref_allele = match.group(2)
            alt_allele = match.group(3)
            flip_sign = 1
            if effect_allele is None:
                effect_allele = match.group(4)
            elif effect_allele != match.group(4):
                flip_sign = -1
            if match.group(7) != "NA":
                sample_size = int(match.group(5))
                if model_name == "glm.logistic.hybrid":
                    betas.append(flip_sign * math.log(float(match.group(6))))
                else:
                    betas.append(flip_sign * float(match.group(6)))
                variances.append(pow(float(match.group(7)), 2))
        else:
            raise ValueError(
                "plink2 regression result does not match expected format: {}".format(
                    line
                )
            )
    return format_mi_result(
        varid,
        effect_allele,
        ref_allele if ref_allele != effect_allele else alt_allele,
        sample_size,
        betas,
        variances,
        len(lines),
    )


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
        out.writelines(
            "varid\teffect_allele\tother_allele\teffect\twithin_variance\tbetween_variance\ttotal_variance\tdof\ttstat\n"
        )
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
