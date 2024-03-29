import math
import os
import pathlib
import runpy

import numpy as np
import pandas as pd
import pytest
import snakemake.script as sms
from pandas.testing import assert_frame_equal
from snakemake.io import Namedlist


@pytest.fixture
def common_tmpdir(tmp_path):
    """
    Create a single temporary output
    directory for all fixtures
    """
    return tmp_path


@pytest.fixture
def variant_count():
    """
    Fix variant count for all test data
    """
    return 5


@pytest.fixture
def sample_size():
    """
    Fix reported sample size for all test data
    """
    return 100


@pytest.fixture
def chromosomes(variant_count):
    """
    Get fixed chromosome entries for all test frames
    """
    return ["chr1" for i in range(variant_count)]


@pytest.fixture
def positions(variant_count):
    """
    Get fixed position entries for all test frames
    """
    return [i + 1 for i in range(variant_count)]


@pytest.fixture
def refs(variant_count):
    """
    Get fixed reference alleles for all test frames
    """
    return ["A" for i in range(variant_count)]


@pytest.fixture
def alts(variant_count):
    """
    Get fixed alternate alleles for all test frames
    """
    return ["C" for i in range(variant_count)]


@pytest.fixture
def variant_ids(chromosomes, positions, refs, alts):
    """
    Get fixed variant ID entries for all test frames
    """
    return [
        ":".join([chromosome, str(position), ref, alt])
        for (chromosome, position, ref, alt) in zip(chromosomes, positions, refs, alts)
    ]


@pytest.fixture
def first_input_a1(refs):
    """
    Shuffle order of reported effect allele for first
    comparison file
    """
    return refs


@pytest.fixture
def second_input_a1(refs, alts):
    """
    Shuffle order of reported effect allele for second
    comparison file
    """
    return [(refs[0] if i % 2 == 0 else alts[0]) for i in range(len(refs))]


@pytest.fixture
def first_input_df_plink2_linear(
    chromosomes, positions, variant_ids, refs, alts, first_input_a1, sample_size
):
    """
    Create a simple DataFrame with
    plink2 linear style output contents

    Observations:
    - many columns are ignored: A1_CT, A1_FREQ, TEST, P
    - A1 needs to be different between test files to test allele matching
    - need test of invariant betas
    - need test of NA entry handling, detected as SE == np.nan
    """
    df = pd.DataFrame(
        {
            "#CHROM": chromosomes,
            "POS": positions,
            "ID": variant_ids,
            "REF": refs,
            "ALT": alts,
            "A1": first_input_a1,
            "A1_CT": "NA",
            "A1_FREQ": "NA",
            "TEST": "ADD",
            "OBS_CT": sample_size,
            "BETA": [1.01, 0.55, -0.33, 0.33, 0.01],
            "SE": ["0.45", "0.44", "NA", "0.22", "NA"],
            "P": "NA",
        }
    )
    return df


@pytest.fixture
def first_input_df_plink2_logistic(first_input_df_plink2_linear):
    """
    Create an input plink2 glm.logistic.hybrid file by slightly
    modifying the first linear regression test input file
    """
    lindf = first_input_df_plink2_linear
    df = pd.DataFrame(
        {
            "#CHROM": lindf["#CHROM"],
            "POS": lindf["POS"],
            "ID": lindf["ID"],
            "REF": lindf["REF"],
            "ALT": lindf["ALT"],
            "A1": lindf["A1"],
            "A1_CT": lindf["A1_CT"],
            "A1_CASE_CT": lindf["A1_CT"],
            "A1_CTRL_CT": lindf["A1_CT"],
            "A1_FREQ": lindf["A1_FREQ"],
            "A1_CASE_FREQ": lindf["A1_FREQ"],
            "A1_CTRL_FREQ": lindf["A1_FREQ"],
            "TEST": lindf["TEST"],
            "OBS_CT": lindf["OBS_CT"],
            "OR": [math.exp(beta) for beta in lindf["BETA"]],
            "LOG(OR)_SE": lindf["SE"],
            "P": lindf["P"],
        }
    )
    return df


@pytest.fixture
def second_input_df_plink2_linear(
    chromosomes, positions, variant_ids, refs, alts, second_input_a1, sample_size
):
    """
    Create a simple DataFrame with
    plink2 linear style output contents

    Observations:
    - many columns are ignored: A1_CT, A1_FREQ, TEST, P
    - A1 needs to be different between test files to test allele matching
    - need test of invariant betas
    - need test of NA entry handling, detected as SE == np.nan
    """
    df = pd.DataFrame(
        {
            "#CHROM": chromosomes,
            "POS": positions,
            "ID": variant_ids,
            "REF": refs,
            "ALT": alts,
            "A1": second_input_a1,
            "A1_CT": "NA",
            "A1_FREQ": "NA",
            "TEST": "ADD",
            "OBS_CT": sample_size,
            "BETA": [0.55, -0.45, -0.36, -0.36, 0.06],
            "SE": ["0.35", "0.34", "0.21", "NA", "NA"],
            "P": "NA",
        }
    )
    return df


@pytest.fixture
def second_input_df_plink2_logistic(second_input_df_plink2_linear):
    """
    Create an input plink2 glm.logistic.hybrid file by slightly
    modifying the second linear regression test input file
    """
    lindf = second_input_df_plink2_linear
    df = pd.DataFrame(
        {
            "#CHROM": lindf["#CHROM"],
            "POS": lindf["POS"],
            "ID": lindf["ID"],
            "REF": lindf["REF"],
            "ALT": lindf["ALT"],
            "A1": lindf["A1"],
            "A1_CT": lindf["A1_CT"],
            "A1_CASE_CT": lindf["A1_CT"],
            "A1_CTRL_CT": lindf["A1_CT"],
            "A1_FREQ": lindf["A1_FREQ"],
            "A1_CASE_FREQ": lindf["A1_FREQ"],
            "A1_CTRL_FREQ": lindf["A1_FREQ"],
            "TEST": lindf["TEST"],
            "OBS_CT": lindf["OBS_CT"],
            "OR": [math.exp(beta) for beta in lindf["BETA"]],
            "LOG(OR)_SE": lindf["SE"],
            "P": lindf["P"],
        }
    )
    return df


@pytest.fixture
def malformed_input_df(chromosomes, positions, variant_ids, refs, alts, sample_size):
    """
    Create a simple DataFrame that, for simplicity,
    doesn't match the format of any of the input types
    accepted by combine_analysis_output
    """
    df = pd.DataFrame(
        {
            "#CHROM": chromosomes,
            "POS": positions,
            "ID": variant_ids,
            "REF": refs,
            "ALT": alts,
            "P": "NA",
        }
    )
    return df


@pytest.fixture
def input_file1_plink2_linear(first_input_df_plink2_linear, common_tmpdir):
    """
    Create a temporary filename and record input data in it
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_plink2_linear_f1.tsv")
    first_input_df_plink2_linear.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def truncated_input(first_input_df_plink2_linear, common_tmpdir):
    """
    Create a temporary filename and record in it the same input
    as plink2 linear regression file 1, but with fewer lines to
    hopefully trigger a jagged input file error
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_plink2_linear_f1_jagged.tsv")
    truncated_df = first_input_df_plink2_linear.iloc[
        [i for i in range(len(first_input_df_plink2_linear) - 1)]
    ]
    truncated_df.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def input_file2_plink2_linear(second_input_df_plink2_linear, common_tmpdir):
    """
    Create a temporary filename and record input data in it
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_plink2_linear_f2.tsv")
    second_input_df_plink2_linear.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def input_file1_plink2_logistic(first_input_df_plink2_logistic, common_tmpdir):
    """
    Create a temporary filename and record input data in it
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_plink2_logistic_f1.tsv")
    first_input_df_plink2_logistic.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def input_file2_plink2_logistic(second_input_df_plink2_logistic, common_tmpdir):
    """
    Create a temporary filename and record input data in it
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_plink2_logistic_f2.tsv")
    second_input_df_plink2_logistic.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def malformed_input(malformed_input_df, common_tmpdir):
    """
    Create a temporary filename and record malformed input data in it
    """
    td = common_tmpdir
    fn = str(td / "combine_analysis_output_malformed_input.tsv")
    malformed_input_df.to_csv(fn, sep="\t", index=False)
    return fn


@pytest.fixture
def output_filename(common_tmpdir):
    """
    Get a temporary filename target for run
    """
    td = common_tmpdir
    return str(td / "combine_analysis_output_results.tsv")


@pytest.fixture
def expected_df_plink2_linear(variant_ids, refs, alts):
    """
    for assertion, record expected MI result
    for plink2 linear example
    """
    df = pd.DataFrame(
        {
            "varid": variant_ids,
            "effect_allele": refs,
            "other_allele": alts,
            "effect": [0.780, 0.500, np.nan, np.nan, np.nan],
            "within_variance": [0.1625, 0.1546, np.nan, np.nan, np.nan],
            "between_variance": [0.0529, 0.0025, np.nan, np.nan, np.nan],
            "total_variance": [0.24185, 0.15835, np.nan, np.nan, np.nan],
            "dof": [8.120801, 89.098287, np.nan, np.nan, np.nan],
            "tstat": [1.5860671, 1.256496, np.nan, np.nan, np.nan],
        }
    )
    return df


@pytest.fixture
def expected_df_plink2_logistic(expected_df_plink2_linear):
    """
    for assertion, record expected MI result
    for plink2 logistic example. based on how
    this was constructed, it should literally
    be the same thing as the linear case
    """
    return expected_df_plink2_linear


@pytest.mark.parametrize(
    "input_file1,input_file2,output_file,params_tool,params_model,expected_df",
    [
        (
            pytest.lazy_fixture("input_file1_plink2_linear"),
            pytest.lazy_fixture("input_file2_plink2_linear"),
            pytest.lazy_fixture("output_filename"),
            "plink2",
            "glm.linear",
            pytest.lazy_fixture("expected_df_plink2_linear"),
        ),
        (
            pytest.lazy_fixture("input_file1_plink2_logistic"),
            pytest.lazy_fixture("input_file2_plink2_logistic"),
            pytest.lazy_fixture("output_filename"),
            "plink2",
            "glm.logistic.hybrid",
            pytest.lazy_fixture("expected_df_plink2_logistic"),
        ),
    ],
)
def test_combine_analysis_output(
    input_file1, input_file2, output_file, params_tool, params_model, expected_df
):
    """
    Given a combination of input settings, test that combine_analysis_output
    can generate the expected output
    """
    snakemake_input = Namedlist(fromdict={"0": input_file1, "1": input_file2})
    snakemake_output = Namedlist(fromdict={"0": output_file})
    snakemake_params = Namedlist(fromdict={"tool": params_tool, "model": params_model})
    smk = sms.Snakemake(
        snakemake_input,
        snakemake_output,
        snakemake_params,
        Namedlist(),
        1,
        Namedlist(),
        Namedlist(),
        {},
        "",
        [],
    )
    ## override depickling
    smk.input = snakemake_input
    smk.output = snakemake_output
    smk.params = snakemake_params
    ## evaluate
    pathlib.Path(os.path.dirname(smk.output[0])).mkdir(parents=True, exist_ok=True)
    runpy.run_path(
        "workflow/scripts/combine_analysis_output.py", init_globals={"snakemake": smk}
    )
    ## load observed output data frame
    assert pathlib.Path(smk.output[0]).is_file()
    observed = pd.read_table(smk.output[0], sep="\t")
    assert_frame_equal(observed, expected_df)


def test_script_insulation():
    """
    Expect that embedded script exits gracefully and without
    error if no snakemake object is defined
    """
    try:
        runpy.run_path("workflow/scripts/combine_analysis_output.py")
    except Exception:
        pytest.fail("evaluation of script with no snakemake object failed")


@pytest.fixture
def dummy_smk(input_file1_plink2_linear, input_file2_plink2_linear, output_filename):
    """
    Construct a simple snakemake object for use
    with error condition testers
    """
    snakemake_input = Namedlist(
        fromdict={"0": input_file1_plink2_linear, "1": input_file2_plink2_linear}
    )
    snakemake_output = Namedlist(fromdict={"0": output_filename})
    snakemake_params = Namedlist(fromdict={"tool": "plink2", "model": "glm.linear"})
    smk = sms.Snakemake(
        snakemake_input,
        snakemake_output,
        snakemake_params,
        Namedlist(),
        1,
        Namedlist(),
        Namedlist(),
        {},
        "",
        [],
    )
    ## override depickling
    smk.input = snakemake_input
    smk.output = snakemake_output
    smk.params = snakemake_params
    return smk


def test_invalid_model_handler(dummy_smk):
    """
    Expect that an invalid model specification results
    in a raised exception
    """
    dummy_smk.params = Namedlist(
        fromdict={"model": "quadratic", "tool": dummy_smk.params["tool"]}
    )
    with pytest.raises(ValueError):
        runpy.run_path(
            "workflow/scripts/combine_analysis_output.py",
            init_globals={"snakemake": dummy_smk},
        )


def test_invalid_tool_handler(dummy_smk):
    """
    Expect that an invalid tool specification results
    in a raised exception
    """
    dummy_smk.params = Namedlist(
        fromdict={"model": dummy_smk.params["model"], "tool": "plink42"}
    )
    with pytest.raises(ValueError):
        runpy.run_path(
            "workflow/scripts/combine_analysis_output.py",
            init_globals={"snakemake": dummy_smk},
        )


@pytest.mark.parametrize(
    "toolname,modelname",
    [("plink2", "glm.linear"), ("plink2", "glm.logistic")],
)
def test_malformed_line_handler(dummy_smk, malformed_input, toolname, modelname):
    """
    Expect that an invalid input line results
    in a raised exception for all tools/models
    """
    dummy_smk.input[0] = malformed_input
    dummy_smk.params = Namedlist(fromdict={"tool": toolname, "model": modelname})
    with pytest.raises(ValueError):
        runpy.run_path(
            "workflow/scripts/combine_analysis_output.py",
            init_globals={"snakemake": dummy_smk},
        )


def test_jagged_file_handler(dummy_smk, truncated_input):
    """
    Expect that code can detect when input files don't have
    the same number of lines
    """
    dummy_smk.input[0] = truncated_input
    with pytest.raises(ValueError):
        runpy.run_path(
            "workflow/scripts/combine_analysis_output.py",
            init_globals={"snakemake": dummy_smk},
        )
