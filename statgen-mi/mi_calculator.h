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

#ifndef __MI_CALCULATOR_H__
#define __MI_CALCULATOR_H__

#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream>

namespace MI {
  class mi_calculator {
  public:
    mi_calculator() {}
    ~mi_calculator() throw() {}

    std::string report_result(const std::string &rsid,
			      std::vector<double> &betas,
			      std::vector<double> &stderrs,
			      const std::vector<std::string> &effect_alleles) const;
  private:
  };
}

#endif //__MI_CALCULATOR_H__
