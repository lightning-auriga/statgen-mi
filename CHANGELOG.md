# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0] - 2022-08-19

### Added

- support for [minimac4](https://genome.sph.umich.edu/wiki/Minimac4) imputation `dose.vcf.gz` probabilities
- support for [plink2](https://www.cog-genomics.org/plink/2.0/) `--glm` regression methods
- [pytest](https://docs.pytest.org/en/7.1.x/) coverage of embedded scripts

### Changed

- workflow is now a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline
- cluster support is now based on [snakemake profiles](https://github.com/Snakemake-Profiles)
- generation of drawn datasets is based on a [modded version](https://github.com/cpalmer718/bcftools) of `bcftools +setGT`
- configuration is handled in yaml under `config/`; see [README](README.md) for details
- dependencies managed with [conda](https://docs.conda.io/en/latest/miniconda.html)
  - eventually, hopefully, the `bcftools +setGT` call will be managed with this as well

### Removed

- support for `impute2`
- support for `snptest`
- both of the above are flagged for readdition later on

[//]: # (- Added)
[//]: # (- Changed)
[//]: # (- Deprecated)
[//]: # (- Removed)
[//]: # (- Fixed)
[//]: # (- Security)
