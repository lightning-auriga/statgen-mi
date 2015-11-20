/*
  Copyright 2015 Lightning Auriga

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

#include "statgen-mi/custom_handlers.h"

void MI::populate_imputed_data_handlers(MI::imputed_input_handler_map &target) {
  imputation_manager im;
  /////////////////////IMPUTE2//////////////////////////
  im.clear();
  im.set_filename_generator(impute2_filename_generator);
  im.add_file_handler(impute2_handle_gen_line);
  im.add_file_handler(impute2_handle_sample_line);
  target["impute2"] = im;
  //////////////////////////////////////////////////////
}

void MI::populate_software_data_handlers(MI::software_input_handler_map &target) {
  program_manager pm;
  /////////////////////SNPTEST//////////////////////////
  pm.clear();
  pm.set_filename_generator(snptest_filename_generator);
  pm.set_genotype_handler(snptest_write_gen_line);
  pm.add_file_handler(snptest_write_sample_line);
  pm.set_command_formatter(snptest_format_command);
  pm.set_result_parser(snptest_process_results);
  target["snptest"] = pm;
  ////////////////////PLINK/////////////////////////////
  pm.clear();
  pm.set_filename_generator(plink_filename_generator);
  pm.set_genotype_handler(plink_write_bed_line);
  pm.add_file_handler(plink_write_bim_line);
  pm.add_file_handler(plink_write_fam_line);
  pm.add_file_handler(plink_write_pheno_line);
  pm.set_command_formatter(plink_format_command);
  pm.set_result_parser(plink_process_results);
  pm.set_extra_file_remover(plink_extra_file_remover);
  target["plink"] = pm;
  ///////////////////EMMAXKIN///////////////////////////
  pm.clear();
  pm.set_filename_generator(emmaxkin_filename_generator);
  pm.set_genotype_handler(emmaxkin_write_tped_line);
  pm.add_file_handler(emmaxkin_write_tfam_line);
  pm.add_file_handler(emmaxkin_write_pheno_line);
  if (!MI::parameters::get_parameter("emmaxkin-emmax-covar-names").empty())
    pm.add_file_handler(emmaxkin_write_covar_line);
  pm.set_command_formatter(emmaxkin_format_command);
  pm.set_result_parser(placeholder_process_results);
  pm.set_extra_file_remover(placeholder_extra_file_remover);
  target["emmaxkin"] = pm;
  //////////////////////////////////////////////////////
}



//#include "statgen-mi/custom_handlers_impute2.cc"
//#include "statgen-mi/custom_handlers_snptest.cc"
//#include "statgen-mi/custom_handlers_plink.cc"
//#include "statgen-mi/custom_handlers_emmaxkin.cc"
//#include "statgen-mi/custom_handlers_placeholder.cc"
