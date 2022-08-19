# statgen-mi

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/lightning-auriga/statgen-mi/tree/default.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/lightning-auriga/statgen-mi/tree/default)

[![codecov](https://codecov.io/gh/lightning-auriga/statgen-mi/branch/default/graph/badge.svg?token=YTKRACNUAN)](https://codecov.io/gh/lightning-auriga/statgen-mi)

## Brief Summary

multiple imputation for association studies

## Authors

* Lightning Auriga (@lightning-auriga)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Clone this repository to your local system, into the place where you want to perform the data analysis.
```
    git clone https://github.com/lightning-auriga/statgen-mi.git
```

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `manifest.tsv` to specify your sample setup.

#### Configuration settings

The following settings are available in primary user configuration under `config/config.yaml`:

- `manifest`: location of primary run manifest file; defaults to `config/manifest.tsv`
- `tools`: **configuration options for association tools supported by the workflow**
  - `bcftools`: **configuration options specific to bcftools**
    - `executable`: command to launch bcftools (see note below)
	- `plugin_path`: path to plugin libraries for this version of bcftools
	- note: this workflow uses `bcftools +setGT` to randomize genotypes from imputed probabilities.
	  this functionality is not actually present in `+setGT`, but has been modded in;
	  see [this repo](https://github.com/lightning-auriga/bcftools/tree/setGT_randomize) for the code.
	  eventually, this will hopefully get migrated somewhere more useful like conda, but a local
	  build suffices for this early stage. recommend placing the bcftools repo at `../bcftools` relative
	  to this workflow
  - `plink2`: **configuration options specific to plink2 --glm methods**
    - `executable`: command to launch plink2. if using conda, this should remain default `plink2`
	- `maxthreads`: maximum number of threads to deploy in a plink2 task
	- `maxmem`: maximum RAM (in MB) supplied to a plink2 task
	- `mi_draws`: number of multiple imputation simulated sets to run for plink2 MI runs
- `imputed_datasets`: **user-defined sets of imputed data that can be selected for analysis**
  - each tag under `imputed_datasets` should be unique, and can be used to refer to the dataset in the manifest
  - each tag should contain under it:
    - `type`: descriptor of imputed file type. currently the only accepted value is `minimac4`
	- `filename`: full path to and name of imputed data file(s). for minimac4: the dose.vcf.gz file(s). if multiple paths are specified in an array, the files will each be processed in turn and concatenated (in order) after run completion
- `regression_models`: **user-defined sets of phenotypes and covariates that can be selected for analysis**
  - each tag under `regression_models` should be unique, and can be used to refer to the model in the manifest
  - each tag should contain under it:
    - `filename`: full path to and name of plink-format phenotype file containing any relevant variables. other variables can also be present
	- `phenotype`: primary outcome for this regression model, as the corresponding header entry in the phenotype file
	- `covariates`: (optional) list of covariates for this regression model, as the corresponding header entry or entries in the phenotype file
	- `model`: descriptor of association type. currently recognized options are `linear` or `logistic`
	- `vif`: (optional) for tools that support this, primarily plink: variance inflation factor cap above which a model is suppressed
- `queue`: **user-defined configuration data for compute queue**
  - `small_partition`: slurm partition (or equivalent for other cluster profiles) for jobs with the following restrictions:
    - max RAM will never exceed 3500M
	- max time will be less than 10 minutes
    - this is exposed to save money, but can just be set to the same value as the below partition if desired
  - `large_partition`: slurm partition (or equivalent for other cluster profiles) for jobs using maximum per-tool analysis settings, with the following additional restrictions:
    - at least 8000M RAM should be available for a task
	- jobs on the partition should be permitted to run at least four hours before being killed

#### Run manifest

Each desired MI run should be configured in a row of the run manifest, by default at `config/manifest.tsv`. The following entries are required for each run:

- `analysis`: unique identifier for this particular run
- `imputed_dataset`: tag for desired imputed dataset to use, as enumerated in `config/config.yaml`
- `tool`: supported association tool for analysis
- `regression_model`: tag for desired regression model to use, as enumerated in `config/config.yaml`


### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --profile /path/to/slurm-profile --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

More information will be added here as it becomes available. For now, draft results are populated under `results/{analysis}/{imputed_dataset}/{tool}/{regression_model}`

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push


## Testing

Testing for embedded snakemake python scripts is in `workflow/scripts/tests` and handled with `pytest`. `snakemake_unit_tests` integration TBD.

## Version History

- see [the changelog](CHANGELOG.md) for details.
