/*
  Copyright 2022 Cameron Palmer

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

#include <string>
#include <iostream>
#include <stdexcept>


//#include "fileinterface/config.h"
#include "fileinterface/fileinterface.h"
#include "statgen-mi/annotations.h"
#include "statgen-mi/boost_random.h"
#include "statgen-mi/cargs.h"
#include "statgen-mi/config_interpreter.h"
#include "statgen-mi/custom_handlers.h"
#include "statgen-mi/mi_calculator.h"
#include "statgen-mi/parameters.h"
#include "statgen-mi/prob_vector.h"
#include "statgen-mi/program_description.h"
#include "statgen-mi/submission_formatter.h"
#include "statgen-mi/transpose_matrix.h"
#include "statgen-mi/utilities.h"


int main(int argc, char **argv) {
  /*
    WORKFLOW IS AS FOLLOWS:
    
    [CUSTOMIZABLE UNIT 1]
    1) read the probabilistic genotypes into memory, or prepare for buffered transpose
    2) read any supplementary files to memory
    
    [CUSTOMIZABLE UNIT 2]
    3) draw datasets and write to file
    4) write other supplementary files for analysis software
    
    [CUSTOMIZABLE UNIT 3]
    5) get a formatted analysis command
    6) submit the draw jobs to queue somewhere
    
    [CUSTOMIZABLE UNIT 4]
    7) in a separate --mi-clean run, combine the draw runs into a single result file and purge old files
  */
  try {
    MI::boost_random::seed();
    MI::cargs arg_parser(argc, argv);
    arg_parser.parse_args(false);

    if (MI::parameters::get_flag("version")) {
      MI::parameters::print_license(std::cout);
      return 0;
    }
    
    std::string imputed_dataset_config = "", submission_queue_config = "", program_config = "";
    MI::config_interpreter interpreter(arg_parser);
    std::map<std::string, MI::imputed_dataset_description> imputed_dataset_stockpile;
    std::map<std::string, MI::imputed_dataset_description>::const_iterator current_imputed_dataset_description;
    std::map<std::string, MI::program_description> program_stockpile;
    std::map<std::string, MI::program_description>::const_iterator current_program_description;
    try {
      //get the configuration files
      imputed_dataset_config = MI::parameters::get_parameter("mi-imputed-dataset-config");
      submission_queue_config = MI::parameters::get_parameter("mi-submission-queue-config");
      program_config = MI::parameters::get_parameter("mi-program-config");
      //parse the configuration files, adding necessary flags to cargs as needed
      interpreter.interpret_imputed_dataset_config(imputed_dataset_config,
						   imputed_dataset_stockpile);
      interpreter.interpret_program_config(program_config,
					   program_stockpile);
      //arg_parser.parse_args(true);
      arg_parser.print_ignored();
    } catch (const std::domain_error &e) {
      std::cerr << "error: could not parse input configuration files" << std::endl;
      MI::parameters::print_help(std::cerr);
      return 42;
    }
    if (argc < 2 || MI::parameters::get_flag("help")) {
      MI::parameters::print_help(std::cout);
      return 0;
    }
		

    //hack the user parameters
    unsigned N = MI::from_string<unsigned>(MI::parameters::get_parameter("mi-draw-number"));
    unsigned baseline_memory = MI::from_string<unsigned>(MI::parameters::get_parameter("mi-baseline-memory"));

    unsigned job_submission_memory_limit = 0;
    unsigned job_submission_time_limit = 0;
    unsigned current_system_memory_limit = 0;

    std::string line = "";
    
    ////////////////////// stockpiles of custom handlers
    MI::imputed_input_handler_map imputed_dataset_line_readers;
    MI::imputed_input_handler_map::iterator current_imputed_dataset_handler;
    MI::software_input_handler_map software_input_line_writers;
    MI::software_input_handler_map::iterator current_software_input_writer;
    
    //have to actually populate the above
    populate_imputed_data_handlers(imputed_dataset_line_readers);
    populate_software_data_handlers(software_input_line_writers);
    
    //now get the user requested types    
    std::string imputed_dataset_type = MI::parameters::get_parameter("imputed-dataset-type");
    std::string program_run_type = MI::parameters::get_parameter("program-run-type");
    //and confirm they actually exist in the various stockpiles
    if ((current_imputed_dataset_description = imputed_dataset_stockpile.find(imputed_dataset_type))
	== imputed_dataset_stockpile.end()) {
      throw std::domain_error("cannot locate dataset type \"" + imputed_dataset_type + "\" within imputed dataset descriptions");
    }
    if ((current_imputed_dataset_handler = imputed_dataset_line_readers.find(imputed_dataset_type))
	== imputed_dataset_line_readers.end()) {
      throw std::domain_error("cannot locate dataset type \"" + imputed_dataset_type + "\" within imputed dataset line reader handlers");
    }
    if ((current_program_description = program_stockpile.find(program_run_type)) == program_stockpile.end()) {
      throw std::domain_error("cannot locate program type \"" + program_run_type + "\" within program descriptions");
    }
    if ((current_software_input_writer = software_input_line_writers.find(program_run_type)) == software_input_line_writers.end()) {
      throw std::domain_error("cannot locate program type \"" + program_run_type + "\" within software input handlers");
    }


    
    
    if (!MI::parameters::get_flag("mi-clean")) {
      MI::annotations annot;
      MI::prob_vector vec;
      ///////////////////////////FIRST YOU READ DATA FROM THE IMPUTED DATA FILES/////////////////////////////////
      //handle the input data with a matrix transposer; will be dealt with later by-line
      MI::matrix_transposer prob_handler(current_imputed_dataset_handler->second.get_filename_generator()(0));
      //handle the remaining input files
      MI::fileinterface_reader *input = 0;
      try {
	for (unsigned all_file_counter = 1;
	     all_file_counter < current_imputed_dataset_handler->second.file_count();
	     ++all_file_counter) {
	  input = MI::reconcile_reader(current_imputed_dataset_handler->second.get_filename_generator()(all_file_counter));
	  //get a line
	  while (input->getline(line)) {
	    current_imputed_dataset_handler->second.get_file_handler(all_file_counter)(line, vec, annot);
	  }
	  input->close();
	  input = 0;
	}
      } catch (...) {
	if (input) delete input;
	input = 0;
	throw;
      }



      /////////////////////////////////////////THEN YOU WRITE DATA FOR MI RUNS/////////////////////////////////////////////////
      std::vector<MI::fileinterface_writer *> drawn_datasets;
      MI::fileinterface_writer *output = 0;
      try {
	for (unsigned current_draw = 1; current_draw <= N; ++current_draw) {
	  output = MI::reconcile_writer(current_software_input_writer->second.get_filename_generator()(0, current_draw));
	  drawn_datasets.push_back(output);
	  output = 0;
	}
	////////if it's the genotype data, do one thing
	//probably gonna wanna warn the user if they're getting the honor of the world's slowest matrix transpose
	if (current_imputed_dataset_description->second.is_snp_major() !=
	    current_program_description->second.is_snp_major())
	  std::cout << "WARNING: input and output do not have same genotype matrix orientation, transpose is happening" << std::endl;
	unsigned current_line = 0;
	while (prob_handler.get_a_line(current_imputed_dataset_handler->second.get_file_handler(0),
				       vec,
				       annot,
				       current_imputed_dataset_description->second.is_snp_major(),
				       current_program_description->second.is_snp_major())) {
	  for (unsigned current_draw = 1; current_draw <= N; ++current_draw) {
	    //if (!current_software_input_writer->second.get_genotype_handler())
	    //throw std::domain_error("genotype handler is null");
	    drawn_datasets.at(current_draw-1)->writeline(current_software_input_writer->second.get_genotype_handler()(vec, annot, current_line));
	  }
	  ++current_line;
	}
	for (unsigned i = 0; i < N; ++i) {
	  drawn_datasets.at(i)->close();
	  delete drawn_datasets.at(i);
	  drawn_datasets.at(i) = 0;
	}
	
	
	
	
	for (unsigned all_file_counter = 0;
	     all_file_counter < current_software_input_writer->second.file_count();
	     ++all_file_counter) {
	  for (unsigned current_draw = 1; current_draw <= N; ++current_draw) {
	    output = MI::reconcile_writer(current_software_input_writer->second.get_filename_generator()(all_file_counter + 1, current_draw));
	    line = current_software_input_writer->second.get_file_handler(all_file_counter)(annot);
	    output->writeline(line);
	    output->close();
	    delete output;
	    output = 0;
	  }
	}
      } catch (...) {
	if (output) delete output;
	for (std::vector<MI::fileinterface_writer *>::iterator iter = drawn_datasets.begin();
	     iter != drawn_datasets.end(); ++iter) {
	  if (*iter) delete *iter;
	}
	throw;
      }


      ////////////////////////TIME TO LAUNCH SOME JOBS////////////////////////////////////
      std::string command = "";
      MI::submission_formatter subformat;
      for (unsigned i = 1; i <= N; ++i) {
	command = current_software_input_writer->second.get_command_formatter()(current_software_input_writer->second.get_filename_generator(),
										i);
	command = current_program_description->second.program_executable() + " " + command;
	command = subformat.get_submission_command(MI::parameters::get_parameter("mi-queue-type"),
						   command,
						   MI::parameters::get_parameter("mi-output-prefix") + ".draw" + MI::to_string<unsigned>(i) + ".output",//output_filename,
						   MI::parameters::get_parameter("mi-output-prefix") + ".draw" + MI::to_string<unsigned>(i) + ".error",//error_filename
						   MI::parameters::get_parameter("mi-baseline-memory"),
						   MI::parameters::get_parameter("mi-time"));
	//std::cout << "I think I should launch the following command:\n\"" << command << "\"" << std::endl;
	if (system(command.c_str())) {
	  throw std::domain_error("unable to launch command \"" + command + "\"");
	}
      }
    } else {
      MI::mi_calculator calc;
      std::vector<MI::fileinterface_reader *> inputs;
      MI::fileinterface_reader *input = 0;
      MI::fileinterface_writer *output = 0;
      try {
	std::vector<std::string> lines;
	lines.resize(N);
	//for each MI draw
	for (unsigned i = 1; i <= N; ++i) {
	  //make sure the job finished successfully
	  if (!MI::successfully_completed(MI::parameters::get_parameter("mi-output-prefix") + ".draw" + MI::to_string<unsigned>(i) + ".output"))
	    throw std::domain_error("detected job failure for draw " + MI::to_string<unsigned>(i));
	  //open the analysis result file
	  input = MI::reconcile_reader(current_software_input_writer->second.get_filename_generator()(current_software_input_writer->second.file_count()+1, i));
	  inputs.push_back(input);
	  input = 0;
	}
	//open the final results file
	output = MI::reconcile_writer(MI::parameters::get_parameter("mi-output-prefix") + ".combined_results.txt");
	//write the final results header
	output->writeline("RSID EFFECTALLELE NDRAWS BETA VAR-WITHIN VAR-BETWEEN VAR-TOTAL DOF TSTAT");
	std::vector<double> betas, stderrs;
	std::vector<std::string> snp_ids, effect_alleles;
	snp_ids.resize(inputs.size(), "");
	betas.resize(inputs.size(), 0.0);
	stderrs.resize(inputs.size(), 0.0);
	effect_alleles.resize(inputs.size(), "");
	while (true) {
	  bool all_files_complete = true;
	  bool all_files_comment = true;
	  for (unsigned i = 0; i < inputs.size(); ++i) {
	    if (!inputs.at(i)->getline(lines.at(i))) {
	      lines.at(i) = "";
	      continue;
	    } else {
	      all_files_complete = false;
	    }
	    if (current_software_input_writer->second.get_result_parser()(lines.at(i), snp_ids.at(i), betas.at(i), stderrs.at(i), effect_alleles.at(0))) {
	      all_files_comment = false;
	    }
	  }
	  if (!all_files_complete && !all_files_comment)
	    output->writeline(calc.report_result(snp_ids.at(0), betas, stderrs, effect_alleles));
	  if (all_files_complete) break;
	}

	for (std::vector<MI::fileinterface_reader *>::iterator iter = inputs.begin(); iter != inputs.end(); ++iter) {
	  if (*iter) {
	    (*iter)->close();
	    delete *iter;
	    *iter = 0;
	  }
	}
	output->close();
	delete output;
	output = 0;
      } catch (...) {
	if (input) delete input;
	for (std::vector<MI::fileinterface_reader *>::iterator iter = inputs.begin(); iter != inputs.end(); ++iter) {
	  if (*iter) delete *iter;
	}
	if (output) delete output;
	throw;
      }



      ///FINALLY: purge files created by this software
      //for each MI draw
      for (unsigned i = 1; i <= N; ++i) {
	//for each file provided by the custom software handler
	for (unsigned j = 0; j <= current_software_input_writer->second.file_count()+1; ++j) {
	  MI::safely_remove(current_software_input_writer->second.get_filename_generator()(j, i));
	}
	//the submission output files, and potentially the error files if using certain queue types
	MI::safely_remove(MI::parameters::get_parameter("mi-output-prefix") + ".draw" + MI::to_string<unsigned>(i) + ".output");
	MI::safely_remove(MI::parameters::get_parameter("mi-output-prefix") + ".draw" + MI::to_string<unsigned>(i) + ".error", false);
	if (current_software_input_writer->second.get_extra_file_remover()) {
	  MI::safely_remove(current_software_input_writer->second.get_extra_file_remover()(i), false);
	}
      }
    }
    std::cout << "all done" << std::endl;
    return 0;
  } catch (const std::domain_error &e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
}
