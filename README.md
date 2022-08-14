# statgen-mi

## Brief Summary

multiple imputation for association studies

## Overview

This README is an automated stub generated from a `cookiecutter` template.
Documentation below reflects the state of the templated project immediately
after creation and may not reflect the current state of the project after
development updates.

## Requirements

  - g++ >= 8.2.0
  - automake/autoconf
  - make >= 4.2
  - git >= 2.28.0
  - nodejs (for commitizen)
  - pre-commit
  - associated linting tools for C++: cppcheck, clang-format
  - [boost headers](https://www.boost.org)
  - [boost program_options](https://www.boost.org/doc/libs/1_75_0/doc/html/program_options.html)
  - [boost filesystem/system](https://www.boost.org/doc/libs/1_75_0/libs/filesystem/doc/index.htm)
  - [boost iostreams](https://www.boost.org/doc/libs/1_74_0/libs/iostreams/doc/index.html)
  - [yaml-cpp](https://github.com/jbeder/yaml-cpp)

## Build

By default, a build process involving a [conda](https://docs.conda.io/en/latest/) environment is supported.

  - if you wish to use `conda` and it's not currently available, you can install it with the instructions [here](https://docs.conda.io/en/latest/miniconda.html)
  - navigate into your project directory (statgen-mi)
  - create the `conda` environment for installation as follows:
  
     `conda env create -f environment.yaml`
  - activate the conda environment:
  
     `conda activate statgen-mi-env`
  - (one time only per environment) install `commitizen`:
  
     `npm install -g commitizen cz-conventional-changelog`
  - (one time only per environment) install `pre-commit` linters:
  
     `pre-commit install`

  - update (create) the necessary `configure` scripts with `autoreconf`:
  
     `autoreconf --force --install`
	 
     - note that this can also be run with `./generate.bash` inside the repo
  - run `configure`:
  
	 `./configure --with-boost=/path/to/miniconda3/envs/statgen-mi-env --with-boost-libdir=/path/to/miniconda3/envs/statgen-mi-env/lib`

	 - if you are planning on installing software to a local directory, run instead `./configure --prefix=/install/dir [...]`
	 - periodically there are some incompatibility issues between `configure` and `conda`. if so, you may need to override
	   some default locations detected by `configure`. for example, you might override the detected compiler with:
	   `CC=gcc CXX=g++ ./configure [...]`
  - run `make CPPFLAGS=""`
	 - this is a non-standard `make` invocation. the reason this is included is because the project
	   is configured to specifically use a `boost` installation in the accompanying `conda` environment.
	   if you'd rather remove `boost` from the conda environment, or ignore it in favor of a system-wide
	   `boost` installation, you can adjust the appropriate `configure` parameters accordingly
	   and instead invoke `make` without any further variable overrides
  - run `make check` to run any `TAP/automake` tests, or the placeholder
     - if you run this command without compiling first, you will again need to override `CPPFLAGS`
	   as follows: `make CPPFLAGS="" check`

  - if desired, run `make install`. if permissions issues are reported, see above for reconfiguring with `./configure --prefix`.
     - as above, if you run installation without compiling first, you will again need to override `CPPFLAGS`
	   as follows: `make CPPFLAGS="" check`
  
## Usage

By default, the final compiled program can be run with

`./statgen-mi.out`

## Version History

14 08 2022: project generated from cookiecutter template

15 11 2015:

- major:
  - Added plink support, linear and logistic traits. PLINK is available [here](http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml)
  - updated example.sh to contain a PLINK-style run type
- minor:
  - fileinterface reconcile_writer now correctly handles .bed{,.gz,.bz2} extensions in filenames

15 11 15 (original release):

Many bugs no doubt exist in this BETA release.  Please send
descriptions of bugs, with sample code invocation as appropriate,
to <cdp2130@cumc.columbia.edu>

This software requires the boost C++ headers.
Currently, the minimum boost version required is 1.46.1

zlib and libbz2 are highly recommended but not strictly required. If included
through ./configure, this will enable files processed by statgen-mi
(primarily imputation genotype files and final result file) to be compressed,
as long as the filenames terminate in ".gz" or ".bz2"

Although this software is designed to interface with other statgen software,
for licensing reasons those software packages are not included with
this release. Currently, statgen-mi is configured to operate with
[IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download)
imputed genotype probabilities and
[SNPTEST2](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)
frequentist analysis tests. statgen-mi needs direct access to the analysis software
(SNPTEST2), so update the --mi-program-config file accordingly
(see examples/mi.program.config for an example configuration).

More information to come as this software approaches maturation.
See examples/example.sh for a test run configuration. Also
./statgen-mi.out --help lists many available options and descriptions.
