UPDATE: 15 nov 2015
major:
-Added plink support, linear and logistic traits. PLINK is available at

http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml

-updated example.sh to contain a PLINK-style run type

minor:
-fileinterface reconcile_writer now correctly handles .bed{,.gz,.bz2} extensions
	       			in filenames




ORIGINAL:
To download for use, select the "statgen-mi.tar.gz" file.

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
IMPUTE2 imputed genotype probabilities and SNPTEST2 frequentist analysis
tests. These software packages can be found, at present, at

https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download

and

https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html

respectively. statgen-mi needs direct access to the analysis software
(SNPTEST2), so update the --mi-program-config file accordingly
(see examples/mi.program.config for an example configuration).

More information to come as this software approaches maturation.
See examples/example.sh for a test run configuration. Also
./statgen-mi.out --help lists many available options and descriptions.







Most recent update: 15 Nov 2015, Cameron Palmer <cdp2130@cumc.columbia.edu>