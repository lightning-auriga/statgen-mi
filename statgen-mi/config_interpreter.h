/*
  Copyright 2022 Lightning Auriga

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

#ifndef STATGEN_MI_CONFIG_INTERPRETER_H__
#define STATGEN_MI_CONFIG_INTERPRETER_H__

#include "statgen-mi/utilities.h"
#include "statgen-mi/cargs.h"
#include "statgen-mi/program_description.h"
#include "fileinterface/fileinterface.h"

namespace statgen_mi {
  class config_interpreter {
  public:
    config_interpreter() {throw std::domain_error("config_interpreter must be provided a cargs object");}
    config_interpreter(cargs &arg_parser)
      : _arg_parser(&arg_parser) {}
    config_interpreter(const config_interpreter &obj)
      : _arg_parser(obj._arg_parser) {}
    ~config_interpreter() throw() {}

    void interpret_imputed_dataset_config(const std::string &filename,
					  std::map<std::string, statgen_mi::imputed_dataset_description> &result);
    void interpret_program_config(const std::string &filename,
				  std::map<std::string, statgen_mi::program_description> &result);
  private:
    bool get_imputed_dataset_chunk(statgen_mi::fileinterface_reader *input,
				   std::string &handle,
				   statgen_mi::imputed_dataset_description &desc) const;
    bool get_program_chunk(statgen_mi::fileinterface_reader *input,
			   std::string &handle,
			   statgen_mi::program_description &desc) const;
    cargs *_arg_parser;
  };
}

#endif //STATGEN_MI_CONFIG_INTERPRETER_H__
