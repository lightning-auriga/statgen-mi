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
def phenotype_df():
    """
    Construct DataFrame representation
    of dummy phenotype file used by all tests
    """
    df = pd.DataFrame(
        {
            "FID": 0,
            "IID": ["A", "B", "C", "D"],
            "phenoname": [1, 2, 3, 4],
            "covar1": [5, 6, 7, 8],
        }
    )
    return df


@pytest.fixture
def phenotype_file(phenotype_df, common_tmpdir):
    """
    Create a valid phenotype file for
    the tests to operate on
    """
    td = common_tmpdir
    out_filename = str(td / "construct_trait_file_phenotype_input.tsv")
    phenotype_df.to_csv(out_filename, sep="\t", index=False)
    return out_filename


@pytest.fixture
def snakemake_input(phenotype_file):
    """
    Create input block featuring valid input target
    """
    return Namedlist(
        fromdict={
            "0": phenotype_file,
        }
    )


@pytest.fixture
def output_filename(common_tmpdir):
    """
    create filename target in tmp space for theoretical output file
    """
    td = common_tmpdir
    res = str(td / "construct_trait_file_phenotype_output.tsv")
    return res


@pytest.fixture
def snakemake_output(output_filename):
    """
    Create output block featuring valid output target
    """
    return Namedlist(
        fromdict={
            "0": output_filename,
        }
    )


@pytest.fixture
def snakemake_params():
    """
    Create a snakemake params block with all valid column requests
    """
    return Namedlist(
        fromdict={
            "phenotype": "phenoname",
            "covariates": "covar1",
        }
    )


@pytest.fixture(
    scope="module",
    params=[("phenoname", ["covar1", "covar2"]), ("badname", ["covar1", "badcovar"])],
)
def snakemake_params_invalid_content(request):
    """
    Create a snakemake params block with some permutation
    of requested column content missing
    """
    return Namedlist(
        fromdict={
            "phenotype": request.param[0],
            "covariates": request.param[1],
        }
    )


@pytest.fixture
def snakemake_object(snakemake_input, snakemake_output, snakemake_params):
    """
    Construct an intact snakemake object with valid entries
    for passing test
    """
    res = sms.Snakemake(
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
    res.input = snakemake_input
    res.output = snakemake_output
    res.params = snakemake_params
    return res


def test_construct_trait_file(snakemake_object, phenotype_df):
    """
    Expect passing test on valid input
    """
    # this is formally only necessary when the output target is in a subdirectory
    # that would usually be constructed by snakemake itself
    pathlib.Path(os.path.dirname(snakemake_object.output[0])).mkdir(
        parents=True, exist_ok=True
    )
    runpy.run_path(
        "workflow/scripts/construct_trait_file.py",
        init_globals={"snakemake": snakemake_object},
    )
    assert pathlib.Path(snakemake_object.output[0]).is_file()
    observed = pd.read_table(snakemake_object.output[0], sep="\t")
    expected = phenotype_df
    assert_frame_equal(observed, expected)
