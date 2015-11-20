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

/////////////////////////////////////////////////PLINK////////////////////////////////////////////////////////////
std::string MI::plink_filename_generator(unsigned index, unsigned draw) {
  switch (index) {
  case 0:
    return MI::parameters::get_parameter("plink-bfile-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("plink-bfile-suffix") + ".bed";
  case 1:
    return MI::parameters::get_parameter("plink-bfile-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("plink-bfile-suffix") + ".bim";
  case 2:
    return MI::parameters::get_parameter("plink-bfile-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("plink-bfile-suffix") + ".fam";
  case 3:
    return MI::parameters::get_parameter("plink-out-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("plink-out-suffix") + ".pheno";
  case 4:
    return MI::parameters::get_parameter("plink-out-prefix") + to_string<unsigned>(draw)
      + MI::parameters::get_parameter("plink-out-suffix")
      + (MI::parameters::get_flag("plink-linear") ? ".assoc.linear" : ".assoc.logistic");
  default:
    throw std::domain_error("plink-filename-generator: invalid index \"" + to_string<unsigned>(index) + "\"");
  }
}
std::string MI::plink_write_bed_line(const MI::prob_vector &vec, const MI::annotations &annot, unsigned index) {
  std::ostringstream o;
  double p1 = 0.0, p2 = 0.0;//, p3 = 0.0;
  unsigned res = 0;
  std::string a1 = annot.get("a1", index);
  std::string a2 = annot.get("a2", index);
  for (unsigned i = 0; i < vec.size(); i += 2) {
    p1 = vec.at(i);
    p2 = vec.at(i+1);
    /*p3 = 1.0 - p1 - p2;
    if (p3 > 1.0) p3 = 1.0;
    if (p3 < 0.0) p3 = 0.0;*/
    res = standard_genotype_draw(p1, p2);
    o << (3-res);
  }
  return o.str();
}
std::string MI::plink_write_bim_line(const MI::annotations &annot) {
  std::ostringstream o;
  std::vector<std::string> rsid = annot.get("rsid");
  std::vector<std::string> pos = annot.get("pos");
  std::vector<std::string> gpos = annot.get("gpos");
  std::vector<std::string> chr = annot.get("chr");
  std::vector<std::string> a1 = annot.get("a1");
  std::vector<std::string> a2 = annot.get("a2");

  if (rsid.empty() ||
      a1.empty() ||
      a2.empty())
    throw std::domain_error("plink-write-bim-line: insufficient snp identifying information (rsid/a1/a2)");
  for (unsigned i = 0; i < rsid.size(); ++i) {
    o << (chr.empty() ? std::string("0") : chr.at(i)) << ' ' << rsid.at(i) << ' '
      << (gpos.empty() ? std::string("0") : gpos.at(i)) << ' '
      << (pos.empty() ? to_string<unsigned>(i+1) : pos.at(i)) << ' '
      << a1.at(i) << ' ' << a2.at(i);
    if (rsid.size() && i < (rsid.size()-1))
      o << '\n';
  }
  return o.str();
}
std::string MI::plink_write_fam_line(const MI::annotations &annot) {
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
std::string MI::plink_write_pheno_line(const MI::annotations &annot) {
  std::string requested_covariates = MI::parameters::get_parameter("plink-covar-names");
  while (requested_covariates.find(",") != std::string::npos)
    requested_covariates[requested_covariates.find(",")] = ' ';
  std::string requested_phenotype = MI::parameters::get_parameter("plink-pheno-name");
  std::ostringstream o;
  o << "FID IID " << requested_phenotype << ' ' << requested_covariates << '\n';
  std::vector<std::string> values;
  std::istringstream strm1(requested_covariates);
  std::string current_covariate = "";
  std::vector<std::string> all_covariates;

  while (strm1 >> current_covariate) {
    all_covariates.push_back(current_covariate);
  }
  for (unsigned i = 0; i < annot.get_count_at("fid"); ++i) {
    o << annot.get("fid", i) << ' ' << annot.get("iid", i) << " ";
    o << annot.get(requested_phenotype, i);
    for (std::vector<std::string>::const_iterator iter = all_covariates.begin(); iter != all_covariates.end(); ++iter) {
      o << ' ' << annot.get(*iter, i);
    }
    if (i != annot.get_count_at("fid")-1) o << '\n';
  }
  return o.str();  
}
std::string MI::plink_format_command(std::string(*filename_generator)(unsigned, unsigned),
				     unsigned draw) {
  std::string command = "";
  command += " --bed " + filename_generator(0, draw) + " --bim " + filename_generator(1, draw)
    + " --fam " + filename_generator(2, draw);
  if (MI::parameters::get_flag("plink-linear") &&
      MI::parameters::get_flag("plink-logistic"))
    throw std::domain_error("plink-format-command: either linear or logistic, not both, should be specified");
  if (MI::parameters::get_flag("plink-linear"))
    command += " --linear";
  else if (MI::parameters::get_flag("plink-logistic"))
    command += " --logistic";
  else
    throw std::domain_error("plink-format-command: neither linear nor logistic mode specified");
  
  command += " --pheno " + filename_generator(3, draw) + " --mpheno 1";
  if (!MI::parameters::get_parameter("plink-covar-names").empty()) {
    command += " --covar-name "
    + MI::parameters::get_parameter("plink-covar-names");
  }
  if (!MI::parameters::get_parameter("plink-vif").empty()) {
    command += " --vif " + MI::parameters::get_parameter("plink-vif");
  }

  command += " --out "
    + MI::parameters::get_parameter("plink-out-prefix") + to_string<unsigned>(draw)
    + MI::parameters::get_parameter("plink-out-suffix");
  return command;  
}
bool MI::plink_process_results(const std::string &result_line,
			       std::string &rsid,
			       double &beta,
			       double &stderr,
			       std::string &effect_allele) {
  if (result_line.find("TEST") != std::string::npos) return false;
  std::istringstream strm1(result_line);
  unsigned pos = 0, nmiss = 0;
  std::string chr = "", test = "";
  if (!(strm1 >> chr >> rsid >> pos >> effect_allele >> test >> nmiss >> beta >> stderr))
    throw std::domain_error("plink-process-results: unable to parse line \"" + result_line + "\"");
  if (MI::parameters::get_flag("plink-logistic")) {
    beta = log(beta);
  }
  stderr = beta / stderr;
  return true;
}
std::string MI::plink_extra_file_remover(unsigned draw) {
  std::string res = "";
  res += MI::parameters::get_parameter("plink-out-prefix") + to_string<unsigned>(draw)
    + MI::parameters::get_parameter("plink-out-suffix") + ".log";
  return res;
}
