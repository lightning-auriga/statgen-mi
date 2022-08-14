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

#include "statgen-mi/custom_handlers.h"

///////////////////////////////SNPTEST////////////////////////////////////////////////////////////////////////////

std::string MI::snptest_filename_generator(unsigned index, unsigned draw) {
  //genotype data first
  switch (index) {
  case 0:
    return MI::parameters::get_parameter("snptest-gen-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("snptest-gen-suffix");
  case 1:
    return MI::parameters::get_parameter("snptest-sample-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("snptest-sample-suffix");
  case 2:
    return MI::parameters::get_parameter("snptest-output-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("snptest-output-suffix");
  default:
    throw std::domain_error("snptest_filename_generator: invalid index \"" + to_string<unsigned>(index) + "\"");
  }
}
 
std::string MI::snptest_write_gen_line(const MI::prob_vector &vec, const MI::annotations &annot, unsigned index) {
  //need chr, rsid, pos, a1, a2 from annot
  std::ostringstream o;
  o << annot.get("chr", index) << ' ' << annot.get("rsid", index) << ' '
    << annot.get("pos", index) << ' ' << annot.get("a1", index) << ' '
    << annot.get("a2", index);
  double p1 = 0.0, p2 = 0.0;//, p3 = 0.0;
  unsigned res = 0;
  for (unsigned i = 0; i < vec.size(); i += 2) {
    p1 = vec.at(i);
    p2 = vec.at(i+1);
    /*p3 = 1.0 - p1 - p2;
    if (p3 > 1.0) p3 = 1.0;
    if (p3 < 0.0) p3 = 0.0;*/
    res = standard_genotype_draw(p1, p2);
    if (res == 1) {
      o << " 1 0 0";
    } else if (res == 2) {
      o << " 0 1 0";
    } else {
      o << " 0 0 1";
    }
  }
  return o.str();
}

std::string MI::snptest_write_sample_line(const MI::annotations &annot) {
  std::string requested_covariates = MI::parameters::get_parameter("snptest-covar");
  while (requested_covariates.find(",") != std::string::npos)
    requested_covariates[requested_covariates.find(",")] = ' ';
  std::string requested_phenotype = MI::parameters::get_parameter("snptest-pheno");
  std::ostringstream o;
  o << "ID_1 ID_2 missing " << requested_phenotype << ' ' << requested_covariates << '\n';
  o << "0 0 0";// P";
  std::vector<std::string> values;
  //sigh. have to determine the types
  bool trait_is_01 = true, trait_is_12 = true, trait_is_integral = true;
  if (annot.get(requested_phenotype + "_type", 0).empty() || MI::parameters::get_flag("snptest-autocorrect-traits")) {
    values = annot.get(requested_phenotype);
    if (values.empty()) throw std::domain_error("snptest: phenotype \"" + MI::parameters::get_parameter("snptest-pheno") + "\" not found");

    snptest_test_trait_type(values, trait_is_01, trait_is_12, trait_is_integral);
    if (trait_is_01 || trait_is_12) {
      o << " B";
    } else if (trait_is_12) {
      o << " P";
    }
  } else {
    o << ' ' << annot.get(requested_phenotype + "_type", 0);
  }

  bool covar_is_01 = true, covar_is_12 = true, covar_is_integral = true;
  std::istringstream strm1(requested_covariates);
  std::string current_covariate = "";
  std::vector<std::string> all_covariates;
  while (strm1 >> current_covariate) {
    all_covariates.push_back(current_covariate);
    if (annot.get(current_covariate + "_type", 0).empty() || MI::parameters::get_flag("snptest-autocorrect-traits")) {
      values = annot.get(current_covariate);
      if (values.empty()) throw std::domain_error("snptest: covariate \"" + current_covariate + "\" not found");
      snptest_test_trait_type(values, covar_is_01, covar_is_12, covar_is_integral);
      if (covar_is_integral) {
	o << " D";
      } else {
	o << " C";
      }
    } else {
      o << ' ' << annot.get(current_covariate + "_type", 0);
    }
  }
  o << '\n';
  //now actually write the values
  for (unsigned i = 0; i < annot.get_count_at("fid"); ++i) {
    o << annot.get("fid", i) << ' ' << annot.get("iid", i) << " 0 ";
    if (trait_is_12) {
      o << (!annot.get(requested_phenotype, i).compare("2") ? "1" : (!annot.get(requested_phenotype, i).compare("1") ? "0" : annot.get(requested_phenotype, i)));
    } else {
      o << annot.get(requested_phenotype, i);
    }
    for (std::vector<std::string>::const_iterator iter = all_covariates.begin(); iter != all_covariates.end(); ++iter) {
      o << ' ' << annot.get(*iter, i);
    }
    if (i != annot.get_count_at("fid")-1) o << '\n';
  }
  return o.str();
}

void MI::snptest_test_trait_type(std::vector<std::string> &values, bool &trait_is_01, bool &trait_is_12, bool &trait_is_integral) {
  trait_is_01 = trait_is_12 = trait_is_integral = true;
  for (std::vector<std::string>::iterator iter = values.begin(); iter != values.end(); ++iter) {
    if (iter->compare("-9") || cicompare(*iter, "na")) {
      *iter = "NA";
      continue;
    }
    if (iter->compare("0") && iter->compare("1")) {
      trait_is_01 = false;
    }
    if (iter->compare("1") && iter->compare("2")) {
      trait_is_12 = false;
    }
    if (iter->find(".") != std::string::npos &&
	iter->find("-") != std::string::npos) {
      trait_is_integral = false;
    }
  }
}

std::string MI::snptest_format_command(std::string(*filename_generator)(unsigned, unsigned), unsigned draw) {
  std::string command = "";
  command += " -data " + filename_generator(0, draw) + " " + filename_generator(1, draw)
    + " -frequentist " + MI::parameters::get_parameter("snptest-frequentist") + " -pheno "
    + MI::parameters::get_parameter("snptest-pheno");
  if (!MI::parameters::get_parameter("snptest-covar").empty())
    command += " -cov_all";
  command += " -method " + MI::parameters::get_parameter("snptest-method") + " -o "
    + filename_generator(2, draw);
  return command;
}

bool MI::snptest_process_results(const std::string &result_line,
				 std::string &rsid,
				 double &beta,
				 double &stderr,
				 std::string &effect_allele) {
  if (result_line.find("#") != std::string::npos ||
      result_line.find_first_not_of(" \t") == std::string::npos ||
      result_line.find("alternate_ids") != std::string::npos) return false;
  std::istringstream strm1(result_line);
  std::string catcher = "";
  //std::cout << "provided line \"" << result_line << "\"" << std::endl;
  if (!(strm1 >> catcher >> rsid))
    throw std::domain_error("snptest: insufficient entries in results line \"" + result_line + "\"");
  for (unsigned i = 0; i < 4; ++i) {
    if (!(strm1 >> effect_allele))
      throw std::domain_error("snptest: insufficient entries in results line \"" + result_line + "\"");
  }
  for (unsigned i = 0; i < 37; ++i) strm1 >> catcher;
  if (!(strm1 >> beta >> stderr)) {
    //throw std::domain_error("snptest: insufficient entries in results line
    beta = stderr = 1.0 / 0.0;
  }
  //std::cout << "got values " << rsid << ' ' << effect_allele << ' ' << beta << ' ' << stderr << std::endl;
  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
