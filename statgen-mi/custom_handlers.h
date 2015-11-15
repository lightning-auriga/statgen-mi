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

#ifndef __MI_CUSTOM_HANDLERS_H__
#define __MI_CUSTOM_HANDLERS_H__

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include "statgen-mi/annotations.h"
#include "statgen-mi/genotype_draw_methods.h"
#include "statgen-mi/parameters.h"
#include "statgen-mi/prob_vector.h"
#include "statgen-mi/utilities.h"

namespace MI {
  class imputation_manager {
  public:
    imputation_manager()
      : _filename_generator(0) {}
    imputation_manager(const imputation_manager &obj)
      : _filename_generator(obj._filename_generator),
      _per_file_handlers(obj._per_file_handlers) {}
    ~imputation_manager() throw() {}

    void set_filename_generator(std::string(*ptr)(unsigned)) {_filename_generator = ptr;}
    std::string(*get_filename_generator())(unsigned) const {return _filename_generator;}
    void add_file_handler(void(*ptr)(const std::string &, MI::prob_vector &, MI::annotations &)) {_per_file_handlers.push_back(ptr);}
    void(*get_file_handler(unsigned index))(const std::string &, MI::prob_vector &, MI::annotations &) const {
      if (index >= _per_file_handlers.size()) throw std::domain_error("imputation_manager: invalid handler index");
      return _per_file_handlers.at(index);
    }
    void clear() {
      _filename_generator = 0;
      _per_file_handlers.clear();
    }
    unsigned file_count() const {return _per_file_handlers.size();}
  private:
    std::string(*_filename_generator)(unsigned);
    typename std::vector<void(*)(const std::string &, MI::prob_vector &, MI::annotations &)> _per_file_handlers;
  };

  class program_manager {
  public:
    program_manager()
      : _filename_generator(0),
      _genotype_handler(0),
      _command_formatter(0),
      _result_parser(0),
      _extra_file_remover(0) {}
    program_manager(const program_manager &obj)
      : _filename_generator(obj._filename_generator),
      _genotype_handler(obj._genotype_handler),
      _nongeno_handlers(obj._nongeno_handlers),
      _command_formatter(obj._command_formatter),
      _result_parser(obj._result_parser),
      _extra_file_remover(obj._extra_file_remover) {}
    ~program_manager() throw() {}

    void set_filename_generator(std::string(*ptr)(unsigned, unsigned)) {_filename_generator = ptr;}
    std::string(*get_filename_generator())(unsigned, unsigned) const {return _filename_generator;}
    void set_genotype_handler(std::string(*ptr)(const MI::prob_vector &, const MI::annotations &, unsigned index)) {_genotype_handler = ptr;}
    std::string(*get_genotype_handler())(const MI::prob_vector &, const MI::annotations &, unsigned index) {return _genotype_handler;}
    void add_file_handler(std::string(*ptr)(const MI::annotations &)) {_nongeno_handlers.push_back(ptr);}
    std::string(*get_file_handler(unsigned index))(const MI::annotations &) const {
      if (index >= _nongeno_handlers.size()) throw std::domain_error("program_manager: invalid handler index");
      return _nongeno_handlers.at(index);
    }
    void set_command_formatter(std::string(*ptr)(std::string(*)(unsigned, unsigned), unsigned)) {_command_formatter = ptr;}
    std::string(*get_command_formatter())(std::string(*)(unsigned, unsigned), unsigned) {return _command_formatter;}

    void set_result_parser(bool(*ptr)(const std::string &, std::string &, double &, double &, std::string &)) {_result_parser = ptr;}
    bool (*get_result_parser())(const std::string &, std::string &, double &, double &, std::string &) {return _result_parser;}

    void set_extra_file_remover(std::string(*ptr)(unsigned)) {_extra_file_remover = ptr;}
    std::string (*get_extra_file_remover())(unsigned) {return _extra_file_remover;}
    void clear() {
      _filename_generator = 0;
      _genotype_handler = 0;
      _nongeno_handlers.clear();
      _command_formatter = 0;
      _result_parser = 0;
      _extra_file_remover = 0;
    }
    unsigned file_count() const {return _nongeno_handlers.size();}
  private:
    std::string(*_filename_generator)(unsigned, unsigned);
    std::string(*_genotype_handler)(const MI::prob_vector &, const MI::annotations &, unsigned index);
    typename std::vector<std::string(*)(const MI::annotations &)> _nongeno_handlers;
    std::string(*_command_formatter)(std::string(*)(unsigned, unsigned), unsigned);
    bool(*_result_parser)(const std::string &, std::string &, double &, double &, std::string &);
    std::string(*_extra_file_remover)(unsigned);
  };
  
  typedef std::map<std::string, imputation_manager> imputed_input_handler_map;
  typedef std::map<std::string, program_manager> software_input_handler_map;
  
  void populate_imputed_data_handlers(imputed_input_handler_map &target);
  void populate_software_data_handlers(software_input_handler_map &target);




  /////////////////////////////////////////////////CUSTOM FUNCTIONS//////////////////////////////////////////

  /*
    editing process:
    1) write one function per imputed data input file
    2) write one function that generates the name of an imputed data file based on an index
    3) create an imputation_manager object with populate_imputed_data_handlers with the
    corresponding tag from the appropriate config file

    similar process for software input file writers, program_manager, and populate_software_data_handlers
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  */
  ///////////////////////////IMPUTE2///////////////////////////////
  std::string impute2_filename_generator(unsigned index);
  void impute2_handle_gen_line(const std::string &line, MI::prob_vector &vec, MI::annotations &annot);
  void impute2_handle_sample_line(const std::string &line, MI::prob_vector &vec, MI::annotations &annot);
  ///////////////////////////SNPTEST///////////////////////////////
  std::string snptest_filename_generator(unsigned index, unsigned draw);
  std::string snptest_write_gen_line(const MI::prob_vector &vec, const MI::annotations &annot, unsigned index);
  std::string snptest_write_sample_line(const MI::annotations &annot);
  std::string snptest_format_command(std::string(*filename_generator)(unsigned, unsigned),
				     unsigned draw);
  bool snptest_process_results(const std::string &result_line,
			       std::string &rsid,
			       double &beta,
			       double &stderr,
			       std::string &effect_allele);
  //supplementary function
  void snptest_test_trait_type(std::vector<std::string> &values, bool &trait_is_01, bool &trait_is_12, bool &trait_is_integral);
  ///////////////////////////PLINK/////////////////////////////////
  std::string plink_filename_generator(unsigned index, unsigned draw);
  std::string plink_write_bed_line(const MI::prob_vector &vec, const MI::annotations &annot, unsigned index);
  std::string plink_write_bim_line(const MI::annotations &annot);
  std::string plink_write_fam_line(const MI::annotations &annot);
  std::string plink_write_pheno_line(const MI::annotations &annot);
  std::string plink_format_command(std::string(*filename_generator)(unsigned, unsigned),
				   unsigned draw);
  bool plink_process_results(const std::string &result_line,
			     std::string &rsid,
			     double &beta,
			     double &stderr,
			     std::string &effect_allele);
  std::string plink_extra_file_remover(unsigned draw);
}


#endif //__MI_CUSTOM_HANDLERS_H__
