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

//////////////////////////////////EMMAXKIN//////////////////////////////////////////////
std::string MI::emmaxkin_filename_generator(unsigned index, unsigned draw) {
  switch (index) {
  case 0:
    return MI::parameters::get_parameter("emmaxkin-tfile-prefix")
      + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("emmaxkin-tfile-suffix")
      + ".tped";
  case 1:
    return MI::parameters::get_parameter("emmaxkin-tfile-prefix")
      + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("emmaxkin-tfile-suffix")
      + ".tfam";
  case 2:
    return MI::parameters::get_parameter("emmaxkin-tfile-prefix")
      + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("emmaxkin-tfile-suffix")
      + ".pheno";
  case 3:
    if (!MI::parameters::get_parameter("emmaxkin-covar-names").empty())
      return MI::parameters::get_parameter("emmaxkin-tfile-prefix")
	+ to_string<unsigned>(draw)
	+ MI::parameters::get_parameter("emmaxkin-tfile-suffix")
	+ ".covar";
  default:
    throw std::domain_error("emmaxkin-filename-generator: invalid index \"" + to_string<unsigned>(index) + "\"");
  }
}
std::string MI::emmaxkin_write_tped_line(const MI::prob_vector &vec, const MI::annotations &annot, unsigned index) {
  std::ostringstream o;
  double p1 = 0.0, p2 = 0.0;//, p3 = 0.0;
  unsigned res = 0;
  std::string a1 = annot.get("a1", index);
  std::string a2 = annot.get("a2", index);
  o << (annot.get("chr", index).empty() ? std::string("1") : annot.get("chr", index)) << ' '
    << annot.get("rsid", index) << ' '
    << (annot.get("gpos", index).empty() ? std::string("0") : annot.get("gpos", index)) << ' '
    << (annot.get("pos", index).empty() ? to_string<unsigned>(index) : annot.get("pos", index));
  for (unsigned i = 0; i < vec.size(); i += 2) {
    p1 = vec.at(i);
    p2 = vec.at(i+1);
    /*p3 = 1.0 - p1 - p2;
    if (p3 > 1.0) p3 = 1.0;
    if (p3 < 0.0) p3 = 0.0;*/
    res = standard_genotype_draw(p1, p2);
    if (res == 1) {
      o << " 1 1"; 
    } else if (res == 2) {
      o << " 1 2";
    } else {
      o << " 2 2";
    }
    //o << (3-res);
  }
  return o.str();
}
std::string MI::emmaxkin_write_tfam_line(const MI::annotations &annot) {
  std::ostringstream o;
  std::vector<std::string> fid = annot.get("fid");
  std::vector<std::string> iid = annot.get("iid");
  std::vector<std::string> pat = annot.get("pat");
  std::vector<std::string> mat = annot.get("mat");
  std::vector<std::string> sex = annot.get("sex");
  //std::vector<std::string> pheno = annot.get("pheno");
  if (iid.empty())
    throw std::domain_error("plink-write-fam-line: no individual IDs present");
  for (unsigned i = 0; i < iid.size(); ++i) {
    o << (fid.empty() ? iid.at(i) : fid.at(i)) << ' ' << iid.at(i) << ' '
      << (pat.empty() ? std::string("0") : pat.at(i)) << ' '
      << (mat.empty() ? std::string("0") : mat.at(i)) << ' '
      << (sex.empty() ? std::string("2") : sex.at(i)) << " -9";
    if (iid.size() && i < (iid.size()-1))
      o << '\n';
  }
  return o.str();
}
std::string MI::emmaxkin_write_pheno_line(const MI::annotations &annot) {
  std::string requested_phenotype = MI::parameters::get_parameter("emmaxkin-emmax-pheno-name");
  std::ostringstream o;
  for (unsigned i = 0; i < annot.get_count_at("fid"); ++i) {
    o << annot.get("fid", i) << ' ' << annot.get("iid", i) << " ";
    o << (annot.get(requested_phenotype, i).compare("-9") ? annot.get(requested_phenotype, i) : std::string("NA"));
    if (i != annot.get_count_at("fid")-1) o << '\n';
  }
  return o.str();  
}

std::string MI::emmaxkin_write_covar_line(const MI::annotations &annot) {
  std::string requested_covariates = MI::parameters::get_parameter("emmaxkin-emmax-covar-names");
  while (requested_covariates.find(",") != std::string::npos)
    requested_covariates[requested_covariates.find(",")] = ' ';
  std::ostringstream o;
  std::istringstream strm1(requested_covariates);
  std::string current_covariate = "";
  std::vector<std::string> all_covariates;

  while (strm1 >> current_covariate) {
    all_covariates.push_back(current_covariate);
  }
  for (unsigned i = 0; i < annot.get_count_at("fid"); ++i) {
    o << annot.get("fid", i) << ' ' << annot.get("iid", i) << "  1 ";
    for (std::vector<std::string>::const_iterator iter = all_covariates.begin(); iter != all_covariates.end(); ++iter) {
      o << ' ' << (annot.get(*iter, i).compare("-9") ? annot.get(*iter, i) : std::string("NA"));
    }
    if (i != annot.get_count_at("fid")-1) o << '\n';
  }
  return o.str();  
}

std::string MI::emmaxkin_format_command(std::string(*filename_generator)(unsigned, unsigned),
					unsigned draw) {
  std::string command = " -v -h ";
  if (MI::parameters::get_flag("emmaxkin-ibs") &&
      MI::parameters::get_flag("emmaxkin-bn")) {
    throw std::domain_error("emmaxkin-format-command: must select either IBS or BN mode");
  } else if (MI::parameters::get_flag("emmaxkin-ibs")) {
    command += "-s ";
  } else if (!MI::parameters::get_flag("emmaxkin-bn")) {
    throw std::domain_error("emmaxkin-format-command: neither IBS nor BN mode selected");
  }
  command += "-d 10 "
    + MI::parameters::get_parameter("emmaxkin-tfile-prefix")
    + to_string<unsigned>(draw)
    + MI::parameters::get_parameter("emmaxkin-tfile-suffix");
  return command;
}
