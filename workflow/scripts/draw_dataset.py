import gzip
import random
import re


def draw_genotype(data_tuple):
    """
    given input triplet probabilities (in fact, just the first
    two of the triplet), draw a random consistent genotype
    """
    p0 = float(data_tuple[0])
    p1 = float(data_tuple[1])
    draw = random.uniform(0, 1)
    res = "1/1"
    if draw < p0:
        res = "0/0"
    elif draw < (p0 + p1):
        res = "0/1"
    return res


def run_draw_dataset(snakemake):
    """
    given an input probabilistic dataset, draw hard call
    output according to the input probability distributions
    """
    random.seed()
    input_filename = snakemake.input[0]
    input_format = snakemake.params["input_format"]
    output_filename = snakemake.output[0]
    output_format = snakemake.params["output_format"]

    # for the moment, just touch output
    if input_format == "minimac4":
        leading_content_pattern = re.compile(
            "([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t"
        )
        gp_content_pattern = re.compile(":([0-9\\.]+),([0-9\\.]+),[0-9\\.]+")
        with gzip.open(input_filename, "rt") as fin, open(output_filename, "w") as fout:
            for input_line in fin:
                if input_line.startswith("#"):
                    fout.write(input_line)
                else:
                    leading_content_match = leading_content_pattern.match(input_line)
                    gp_content_match = gp_content_pattern.findall(input_line)
                    if leading_content_match is None:
                        raise ValueError(
                            "vcf leading content pattern failed to understand input line: {}".format(
                                input_line
                            )
                        )
                    if len(gp_content_match) == 0:
                        raise ValueError(
                            "genotype probability pattern failed to understand input line: {}".format(
                                input_line
                            )
                        )
                    if output_format == "vcf":
                        output_line = "{}\tGT\t{}\n".format(
                            "\t".join(
                                [leading_content_match.group(i + 1) for i in range(8)]
                            ),
                            "\t".join(
                                [
                                    draw_genotype(gp_content_match[i])
                                    for i in range(len(gp_content_match))
                                ]
                            ),
                        )
                    else:
                        raise ValueError(
                            "draw_dataset does not currently support output format {}".format(
                                output_format
                            )
                        )
                    fout.write(output_line)
    else:
        raise ValueError(
            "draw_dataset does not currently support input format {}".format(
                input_format
            )
        )


try:
    snakemake
except NameError:
    pass
else:
    run_draw_dataset(snakemake)  # noqa: F821
