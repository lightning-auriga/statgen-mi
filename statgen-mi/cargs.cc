/*
  Copyright 2015 Cameron Palmer

  This file is part of statgen-mi.
  
  statgen-mi is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  statgen-mi is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with statgen-mi.  If not, see <http://www.gnu.org/licenses/>.
*/

/*!
* \file cargs.cpp
* \brief contains argument parsing class data
*/

#include "statgen-mi/cargs.h"

void MI::cargs::parse_args(bool warn_unused) throw (std::domain_error) {
  //check for each accepted flag

  //set_flag("license", "l", "Print the license to this program", false);
  set_flag("version", "v", "Print the program version and license", false);
  set_flag("help", "h", "Print this help information", false);
	
  set_parameter("mi-draw-number", 
		"midn", 
		"Number of rounds of MI draws to conduct", 
		"10");
  set_flag("mi-clean", 
	   "mic", 
	   "Combine results of MI run into output files", 
	   false);
  set_parameter("mi-baseline-memory",
		"mibm",
		"Starting memory (in GB) of submitted jobs",
		"4");
  set_parameter("mi-time",
		"mit",
		"Time (in hours) to request for a single MI job",
		"127");
  set_parameter("mi-queue-type",
		"mqt",
		"Version of queue (qsub, bsub, none) in -msqc that should be used",
		"qsub");
  set_parameter("mi-output-prefix",
		"miout",
		"Stem for all output files written by the program",
		"mi-results-");
  set_parameter("mi-program-config",
		"mpcc",
		"File listing all permissible MI run types",
		"examples/mi.program.config");
  set_parameter("mi-submission-queue-config",
		"msqc",
		"File listing all available job submission queues",
		"examples/mi.queue.config");
  set_parameter("mi-imputed-dataset-config",
		"midc",
		"File listing all available imputed dataset formats",
		"examples/mi.imputeddataset.config");
  set_parameter("imputed-dataset-type",
		"idt",
		"From -midc, code for which imputed dataset type is in use",
		"impute2");
  set_parameter("program-run-type",
		"prt",
		"From -mpcc, code for which statistical software type is in use",
		"snptest");
  /*
    _configuration.read_file_list(parameters::get_parameter("program-configuration-file"));
    unsigned config_file_counter = 0;
    for (program_configuration::configuration_const_iterator iter = _configuration.begin(); iter != _configuration.end(); ++iter, ++config_file_counter) {
    for (program_configuration::program_const_iterator piter = iter->begin(); piter != iter->end(); ++piter) {
    set_parameter(piter->get_long_name(), piter->get_short_name(), piter->get_description(), "");
    }
    bool all_found = true;
    for (program_configuration::program_const_iterator piter = iter->begin(); all_found && piter != iter->end(); ++piter) {
    if (parameters::get_parameter(piter->get_long_name()).empty()) all_found = false;
    }
    if (all_found) {
    set_flag(_configuration.get_name(config_file_counter) + "-valid-set",
    "",
    "Full parameter set found for group \"" + _configuration.get_name(config_file_counter) + "\"",
    true);
    std::cout << "Detected parameter combination for \"" << _configuration.get_name(config_file_counter) << "\"" << std::endl;
    }
    }
  */
  /*
    set_parameter("1000G-haplotypes-haplo-prefix",
    "1hhp",
    "Prefix to the reference haplotype data files");
    set_parameter("1000G-haplotypes-haplo-suffix",
    "1hhs",
    "Suffix to the reference haplotype data files");
    set_parameter("1000G-haplotypes-legend-prefix",
    "1hlp",
    "Prefix to the reference haplotype legend files");
    set_parameter("1000G-haplotypes-legend-suffix",
    "1hls",
    "Suffix to the reference haplotype legend files");
    set_flag("snp-major",
    "sm",
    "Whether the reference haplotypes are in SNP-major format",
    false);
    set_parameter("chr-first", 
    "cf", 
    "First chromosome to be loaded from file");
    set_parameter("chr-last", 
    "cl", 
    "Last chromosome to be loaded from file");
    set_parameter("number-simulations", 
    "ns", 
    "Number of simulations to run", "100");
    set_parameter("number-simulated-individuals", 
    "nsi", 
    "Number of simulated individuals to draw from reference haplotypes "
    "per simulation");
    set_parameter("snp-inclusion-list", 
    "sil", 
    "List of SNP IDs to be included from input haplotypes");
    set_parameter("snp-exclusion-list", 
    "sel", 
    "List of SNP IDs to be excluded from input haplotypes");
    set_parameter("known-association-results", 
    "kar", 
    "(Headered) list of known GWAS SNPs [rsid a1 a2 beta freq1]");
    set_parameter("geometric-fail-rate", 
    "gfr", 
    "Failure probability of geometric draws", "0.99");
    set_flag("spike-variant", 
    "sv", 
    "Add an interaction to the existing data", 
    false);
    set_parameter("spike-variant-freq1", 
    "svf1", 
    "Frequency of first interacting SNP");
    set_parameter("spike-variant-freq2", 
    "svf2", 
    "Frequency of second interacting SNP");
    set_parameter("spike-variant-variance-explained", 
    "svve", 
    "Variance explained by interaction");
    set_parameter("k", 
    "k", 
    "Number of individuals to be drawn during PAC search", 
    "20");
    set_parameter("T", 
    "T", 
    "Number of times to draw k individuals during PAC search", "1000");
  */
#ifdef HAVE_PTHREAD
  set_parameter("max-threads", 
		"mt", 
		"Number of threads available (when compiled with pthread support)", 
		"2");
#endif //HAVE_PTHREAD
  //open the log file
  //Logfile::openLog(par::_out);
  if (warn_unused)
    print_ignored();
}
void MI::cargs::set_parameter(const std::string &parameter_long, 
			      const std::string &parameter_short, 
			      const std::string &description, 
			      const std::string &default_value) {
  std::string value = default_value;
  if (parameter_long.empty())
    throw std::domain_error("MI::cargs::set_parameter: empty parameter name");
  if (find("--" + parameter_long)) {
    value = get_formatted_value<std::string>("--" + parameter_long);
  } else if (!parameter_short.empty() &&
	     find("-" + parameter_short)) {
    value = get_formatted_value<std::string>("-" + parameter_short);
  }
  parameters::set_parameter(parameter_long, parameter_short, description, value);
}
void MI::cargs::set_flag(const std::string &parameter_long, 
			 const std::string &parameter_short, 
			 const std::string &description, 
			 bool default_value) {
  parameters::set_flag(parameter_long,
		       parameter_short,
		       description,
		       find("--" + parameter_long) || 
		       (!parameter_short.empty() && find("-" + parameter_short)) 
		       ? !default_value 
		       : default_value);
}
void MI::cargs::set_defaults() throw (std::domain_error) {
  _found.resize(num_args(), false);
}
bool MI::cargs::find(const std::string &key) throw (std::domain_error) {
  for (int i = 1; i < num_args(); ++i) {
    if (!std::string(arg_array()[i]).compare(key)) {
      if (_found.at(i)) {
	//throw std::domain_error("User: double use of flag \"" + key + "\"");
      }
      return (_found.at(i) = true);
    }
  }
  return false;
}
void MI::cargs::print_ignored() const throw () {
  for (int i = 1; i < num_args(); ++i) {
    if (!_found.at(i)) {
      std::cout << "unused command line entry: \"" 
		<< std::string(arg_array()[i]) << "\"" << std::endl;
      //Logfile::printLog("Unused command line parameter: \""
      //		+ std::string(arg_array()[i]) + "\"");
    }
  }
}
