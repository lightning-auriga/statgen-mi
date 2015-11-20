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

///////////////////////////////////IMPUTE2////////////////////////////////////////////////////////////////////////
std::string MI::impute2_filename_generator(unsigned index) {
  //genotype data first
  switch (index) {
  case 0:
    return MI::parameters::get_parameter("impute2-gen-filename");
  case 1:
    return MI::parameters::get_parameter("impute2-sample-filename");
  default:
    throw std::domain_error("impute2_filename_generator: invalid index \"" + to_string<unsigned>(index) + "\"");
  }
}

void MI::impute2_handle_gen_line(const std::string &line, MI::prob_vector &vec, MI::annotations &annots) {
  std::string chr = "", rsid = "", pos = "", a1 = "", a2 = "";
  double p1 = 0.0, p2 = 0.0, p3 = 0.0;
  std::istringstream strm1(line);
  if (!(strm1 >> chr >> rsid >> pos >> a1 >> a2))
    throw std::domain_error("invalid impute2 gen line detected: \"" + line + "\"");
  while (strm1 >> p1 >> p2 >> p3) {
    vec.push_back(p1);
    vec.push_back(p2);
  }
  annots.add("chr", chr);
  annots.add("rsid", rsid);
  annots.add("pos", pos);
  annots.add("a1", a1);
  annots.add("a2", a2);
}

void MI::impute2_handle_sample_line(const std::string &line, MI::prob_vector &vec, MI::annotations &annot) {
  //going to have to assume that the sample phenotype/covariate data were preloaded
  //so just ignore the prob_vector
  vec.clear();
  //annot.clear();
  std::string fid = "", iid = "", value = "", missingness = "";
  //std::vector<std::string> pheno_names;
  std::istringstream strm1(line);
  //if it's a header line
  if (line.find("ID_1") != std::string::npos) {//first header line
    if (!(strm1 >> fid >> iid >> missingness))
      throw std::domain_error("impute2_handle_sample_line: cannot parse first header line \"" + line + "\"");
    unsigned val_counter = 1;
    while (strm1 >> value) {
      //std::cout << "logging pair \"header" << val_counter << "\" \"" << value << "\"" << std::endl;
      annot.add("header" + to_string<unsigned>(val_counter), value);
      //pheno_names.push_back(value);
      ++val_counter;
    }
  } else if (line.find("0 0 0") != std::string::npos) {//second header line
    if (!(strm1 >> fid >> iid >> missingness))
      throw std::domain_error("impute2_handle_sample_line: cannot parse second header line \"" + line + "\"");
    unsigned val_counter = 1;
    while (strm1 >> value) {
      //std::cout << "logging type \"" << annot.get("header" + to_string<unsigned>(val_counter),0) << "_type\" \"" << value << "\"" << std::endl;
      annot.add(annot.get("header" + to_string<unsigned>(val_counter), 0) + "_type", value);
      ++val_counter;
    }
  } else {//data line
    if (!(strm1 >> fid >> iid >> missingness))
      throw std::domain_error("impute2_handle_sample_line: cannot parse data line \"" + line + "\"");
    annot.add("fid", fid);
    annot.add("iid", iid);
    annot.add("missingness", missingness);
    unsigned val_counter = 1;
    while (strm1 >> value) {
      //std::cout << "logging phenotype data under value \"" + annot.get("header" + to_string<unsigned>(val_counter), 0) << "\"" << std::endl;
      annot.add(annot.get("header" + to_string<unsigned>(val_counter), 0), value);
      ++val_counter;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
