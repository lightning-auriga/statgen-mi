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

#include "statgen-mi/mi_calculator.h"

std::string statgen_mi::mi_calculator::report_result(const std::string &rsid,
					     std::vector<double> &betas,
					     std::vector<double> &stderrs,
					     const std::vector<std::string> &effect_alleles) const {
  if (betas.size() != stderrs.size() ||
      effect_alleles.size() != betas.size())
    throw std::domain_error("mi_calculator: vector length error");
  unsigned valid_count = 0;
  double mean_beta = 0.0;
  double mean_variance = 0.0;
  for (unsigned i = 0; i < betas.size(); ++i) {
    if (betas.at(i) == betas.at(i)) {
      //std::cout << betas.at(i) << ' ' << stderrs.at(i) << std::endl;
      mean_beta += betas.at(i);
      mean_variance += pow(stderrs.at(i), 2);
      ++valid_count;
    }
  }
  mean_beta /= static_cast<double>(valid_count);
  mean_variance /= static_cast<double>(valid_count);
  double between_variance = 0.0;
  for (unsigned i = 0; i < betas.size(); ++i) {
    if (betas.at(i) == betas.at(i)) {
      between_variance += pow(mean_beta - betas.at(i), 2);
    }
  }
  between_variance /= static_cast<double>(valid_count - 1);

  double total_variance = mean_variance + static_cast<double>(valid_count + 1) / static_cast<double>(valid_count) * between_variance;

  double dof = static_cast<double>(valid_count - 1) * pow(1.0 + 1.0 / (valid_count + 1.0) * mean_variance / between_variance, 2);

  double tstat = mean_beta / sqrt(total_variance);

  std::ostringstream o;
  if (!(o << rsid << ' ' << effect_alleles.at(0) << ' '
	<< valid_count << ' ' << mean_beta << ' '
	<< mean_variance << ' ' << between_variance << ' '
	<< total_variance << ' ' << dof << ' '
	<< tstat))
    throw std::domain_error("mi_calculator: could not format output line");
  return o.str();
}
