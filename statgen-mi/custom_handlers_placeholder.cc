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

//////////////////////////////////UTILITIES/////////////////////////////////////////////
bool MI::placeholder_process_results(const std::string &result_line,
				     std::string &rsid,
				     double &beta,
				     double &stderr,
				     std::string &effect_allele) {
  throw std::domain_error("This program type is not designed to be run in --mi-clean mode");
}
std::string MI::placeholder_extra_file_remover(unsigned draw) {
  throw std::domain_error("This program type is not designed to be run in --mi-clean mode");
}
