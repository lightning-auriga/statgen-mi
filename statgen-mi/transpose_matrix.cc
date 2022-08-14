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

#include "statgen-mi/transpose_matrix.h"


void MI::matrix_transposer::clear_internals() {
  for (subunit_vector::iterator iter = _matrix_subunits.begin();
       iter != _matrix_subunits.end(); ++iter) {
    if (iter->second.get()) {
      iter->second->close();
      remove(iter->first.c_str());
    }
  }
  _matrix_subunits.clear();
  if (_input) {
    _input->close();
    delete _input;
  }
  _per_line_annotations.clear();
}

bool MI::matrix_transposer::transpose_a_matrix(void (*input_converter)(const std::string &, MI::prob_vector &, MI::annotations &),
					       prob_vector &result,
					       unsigned ram_limit,
					       unsigned disk_limit) {
  unsigned max_vector_count = 0;
  std::string line = "";
  prob_vector vec;
  MI::annotations annot;
  //here's the deal:
  /*
    open the file
    find out how many bloody entries there are by running input_converter on a single line
    load lines until you hit ram_limit. if you hit it, dump to file and start file transpose mode
    otherwise, just keep the thing sitting in memory and give it to the caller in lines
   */
  /////////////////we need a single line from the data. Did we already do the transposition?
  if (_matrix_subunits.empty()) {
    //were we able to just load the bloody thing into RAM?
    if (!_single_matrix.empty()) {
      return get_a_line_from_single_matrix(result);
    } else { //we have to actually deal with the data this time
      while (_input->getline(line)) {
	(*input_converter)(line, vec, annot);
	_single_matrix.push_back(vec);
	_per_line_annotations.push_back(annot);
	max_vector_count = static_cast<unsigned>(ceil(static_cast<double>(ram_limit * 1000000)/(2.5*static_cast<double>(vec.size())))) * 1000;
	//check the total size of the in-memory matrix
	if (_single_matrix.max_size() >= max_vector_count &&
	    _single_matrix.size() == _single_matrix.max_size()) {
	  write_current_matrix();
	}
      }
    }
  }
  //now we know the data are loaded
  if (!_single_matrix.empty()) {//it's in a single in-RAM matrix
    return get_a_line_from_single_matrix(result);
  } else {//it's in the files
    return get_a_line_from_multiple_files(result);
  }
}

bool MI::matrix_transposer::get_a_line_from_single_matrix(prob_vector &vec) {
  //this one's easier. get the values from a column in a single prob_vector
  //MI::prob_vector vec;
  vec.clear();
  if (_current_column == _single_matrix.at(0).size()) return false;
  for (unsigned i = 0; i < _single_matrix.size(); ++i) {
    vec.push_back(_single_matrix.at(i).at(_current_column));
  }
  ++_current_column;
  //(*output_converter)(vec, result);
  return true;
}

bool MI::matrix_transposer::get_a_line_from_multiple_files(prob_vector &vec) {
  std::string line = "";
  //MI::prob_vector vec;
  double p1 = 0.0, p2 = 0.0;//  std::pair<double, double> val;
  
  for (subunit_vector::iterator iter = _matrix_subunits.begin(); iter != _matrix_subunits.end(); ++iter) {
    //get a line
    if (!iter->second->getline(line)) return false;
    std::istringstream strm1(line);
    while (strm1 >> p1 >> p2) {
      vec.push_back(p1);
      vec.push_back(p2);
    }
  }
  //(*output_converter)(vec, result);
  return true;
}

bool MI::matrix_transposer::get_a_line(void (*input_converter)(const std::string &, MI::prob_vector &, MI::annotations &),
				       prob_vector &result,
				       annotations &result_annotation,
				       bool input_is_snp_major,
				       bool output_is_snp_major) {
  std::string line = "";
  MI::prob_vector vec;
  result.clear();
  //MI::annotations annot;
  if (!_input) _input = MI::reconcile_reader(get_filename());
  if (input_is_snp_major == output_is_snp_major) {
    if (_input->getline(line)) {
      (*input_converter)(line, result, result_annotation);
      _per_line_annotations.push_back(result_annotation);
      //(*output_converter)(vec, result);
      return true;
    } else return false;
  } else {
    result_annotation.clear();
    return transpose_a_matrix(input_converter,
			      result);
  }
}

void MI::matrix_transposer::write_current_matrix() {
  std::string a_filename = "";
  std::ifstream input;
  while (true) {
    a_filename = boost::lexical_cast<std::string>((boost::uuids::random_generator())()) + ".gz";
    input.open(a_filename.c_str());
    if (input.is_open()) {
      input.close();
      input.clear();
    } else {
      input.clear();
      break;
    }
  }
  MI::fileinterface_writer *output = 0;
  try {
    output = MI::reconcile_writer(a_filename.c_str());
    for (unsigned i = 0; i < _single_matrix.at(0).size(); i += 2) {
      std::string line = "";
      for (unsigned j = 0; j < _single_matrix.size(); ++j) {
	if (j) line += " ";
	line += to_string<double>(_single_matrix.at(j).at(i)) + " "
	  + to_string<double>(_single_matrix.at(j).at(i+1));
      }
      output->writeline(line);
    }
    output->close();
    delete output;
    output = 0;
  } catch (...) {
    if (output) {delete output; output = 0;}
    throw;
  }
  _single_matrix.clear();
  _current_column = 0;
  boost::shared_ptr<fileinterface_reader_type> stored_input(new fileinterface_reader_type);
  stored_input->open(a_filename.c_str());
  _matrix_subunits.push_back(std::make_pair(a_filename, stored_input));
}
