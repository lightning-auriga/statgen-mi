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

#include "statgen-mi/config_interpreter.h"

void statgen_mi::config_interpreter::interpret_imputed_dataset_config(const std::string &filename,
							      std::map<std::string, statgen_mi::imputed_dataset_description> &result) {
  statgen_mi::fileinterface_reader *input = 0;
  std::string line = "", handle = "";
  statgen_mi::imputed_dataset_description desc;
  try {
    input = statgen_mi::reconcile_reader(filename);
    while (get_imputed_dataset_chunk(input, handle, desc)) {
      result[handle] = desc;
    }
    input->close();
    input = 0;
  } catch (...) {
    if (input) delete input;
    throw;
  }
}
void statgen_mi::config_interpreter::interpret_program_config(const std::string &filename,
						      std::map<std::string, statgen_mi::program_description> &result) {
  statgen_mi::fileinterface_reader *input = 0;
  std::string line = "", handle = "";
  statgen_mi::program_description desc;
  try {
    input = statgen_mi::reconcile_reader(filename);
    while (get_program_chunk(input, handle, desc)) {
      result[handle] = desc;
    }
    input->close();
    input = 0;
  } catch (...) {
    if (input) delete input;
    throw;
  }
}

bool statgen_mi::config_interpreter::get_imputed_dataset_chunk(statgen_mi::fileinterface_reader *input,
						       std::string &handle,
						       statgen_mi::imputed_dataset_description &desc) const {
  std::string line = "", tag = "", value = "", flag_description = "";
  bool within_fragment = false;
  while (input->getline(line)) {
    line = line.substr(0, line.find("#"));
    if (line.empty() || line.find_first_not_of(" \t,") == std::string::npos) continue;
    //line = line.substr(line.find_first_not_of(" \t,"));
    std::istringstream strm1(line);
    if (!(strm1 >> tag)) continue;
    if (cicompare("begin", tag)) {
      if (within_fragment) throw std::domain_error("config_interpreter::imputed_dataset_config: restarted within existing set");
      within_fragment = true;
      if (!(strm1 >> handle)) throw std::domain_error("config_interpreter::imputed_dataset_config: no name for fragment");
      desc.set_program_handle(handle);
    } else if (cicompare("end", tag)) {
      if (!within_fragment) throw std::domain_error("config_interpreter::imputed_dataset_config: stopped outside of set");
      within_fragment = false;
      return true;
    } else if (cicompare("snp_major", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::imputed_dataset_config: no value for tag snp_major");
      if (cicompare(value, "t") ||
	  cicompare(value, "true") ||
	  cicompare(value, "1") ||
	  cicompare(value, "y") ||
	  cicompare(value, "yes")) {
	desc.set_snp_major(true);
      } else if (cicompare(value, "f") ||
		 cicompare(value, "false") ||
		 cicompare(value, "0") ||
		 cicompare(value, "n") ||
		 cicompare(value, "no")) {
	desc.set_snp_major(false);
      } else {
	throw std::domain_error("config_interpreter::imputed_dataset_config: unrecognized value for tag snp_major: \"" + value + "\"");
      }
    } else if (cicompare("parameter", tag) ||
	       cicompare("param", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::imputed_dataset_config: no value for tag parameter");
      //get the quoted description
      if (line.find("\"") == std::string::npos ||
	  line.find("\"") == line.rfind("\"")) throw std::domain_error("config_interpreter::imputed_dataset_config: parameter line requires quoted description");
      flag_description = line.substr(line.find("\"") + 1, line.rfind("\"") - line.find("\"") - 1);

      _arg_parser->set_parameter(desc.program_handle() + "-" + value,
				 value,
				 flag_description,
				 "");
    } else if (cicompare("flag", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::imputed_dataset_config: no value for tag flag");
      //get the quoted description
      if (line.find("\"") == std::string::npos ||
	  line.find("\"") == line.rfind("\"")) throw std::domain_error("config_interpreter::imputed_dataset_config: flag line requires quoted description");
      flag_description = line.substr(line.find("\"") + 1, line.rfind("\"") - line.find("\"") - 1);

      _arg_parser->set_flag(desc.program_handle() + "-" + value,
			    value,
			    flag_description,
			    false);
    } else {
      throw std::domain_error("config_interpreter::imputed_dataset_config: unrecognized tag: \"" + tag + "\"");
    }
  }
  return false;
}
bool statgen_mi::config_interpreter::get_program_chunk(statgen_mi::fileinterface_reader *input,
					       std::string &handle,
					       statgen_mi::program_description &desc) const {
  std::string line = "", tag = "", value = "", flag_description = "";
  bool within_fragment = false;
  while (input->getline(line)) {
    line = line.substr(0, line.find("#"));
    if (line.empty() || line.find_first_not_of(" \t,") == std::string::npos) continue;
    //line = line.substr(line.find_first_not_of(" \t,"));
    std::istringstream strm1(line);
    if (!(strm1 >> tag)) continue;
    if (cicompare("begin", tag)) {
      if (within_fragment) throw std::domain_error("config_interpreter::program_config: restarted within existing set");
      within_fragment = true;
      if (!(strm1 >> handle)) throw std::domain_error("config_interpreter::program_config: no name for fragment");
      desc.set_program_handle(handle);
    } else if (cicompare("end", tag)) {
      if (!within_fragment) throw std::domain_error("config_interpreter::program_config: stopped outside of set");
      within_fragment = false;
      return true;
    } else if (cicompare("snp_major", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::program_config: no value for tag snp_major");
      if (cicompare(value, "t") ||
	  cicompare(value, "true") ||
	  cicompare(value, "1") ||
	  cicompare(value, "y") ||
	  cicompare(value, "yes")) {
	desc.set_snp_major(true);
      } else if (cicompare(value, "f") ||
		 cicompare(value, "false") ||
		 cicompare(value, "0") ||
		 cicompare(value, "n") ||
		 cicompare(value, "no")) {
	desc.set_snp_major(false);
      } else {
	throw std::domain_error("config_interpreter::program_config: unrecognized value for tag snp_major: \"" + value + "\"");
      }
    } else if (cicompare("parameter", tag) ||
	       cicompare("param", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::program_config: no value for tag parameter");
      //get the quoted description
      if (line.find("\"") == std::string::npos ||
	  line.find("\"") == line.rfind("\"")) throw std::domain_error("config_interpreter::program_config: parameter line requires quoted description");
      flag_description = line.substr(line.find("\"") + 1, line.rfind("\"") - line.find("\"") - 1);

      _arg_parser->set_parameter(desc.program_handle() + "-" + value,
				 value,
				 flag_description,
				 "");
    } else if (cicompare("flag", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::program_config: no value for tag flag");
      //get the quoted description
      if (line.find("\"") == std::string::npos ||
	  line.find("\"") == line.rfind("\"")) throw std::domain_error("config_interpreter::program_config: flag line requires quoted description");
      flag_description = line.substr(line.find("\"") + 1, line.rfind("\"") - line.find("\"") - 1);

      _arg_parser->set_flag(desc.program_handle() + "-" + value,
			    value,
			    flag_description,
			    false);
    } else if (cicompare("executable", tag)) {
      if (!(strm1 >> value)) throw std::domain_error("config_interpreter::program_config: no value for tag executable");
      desc.set_program_executable(value);
    } else {
      throw std::domain_error("config_interpreter::program_config: unrecognized tag: \"" + tag + "\"");
    }
  }
  return false;  
}
