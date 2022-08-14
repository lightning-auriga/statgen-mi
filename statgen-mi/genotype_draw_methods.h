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

#ifndef STATGEN_MI_GENOTYPE_DRAW_METHODS_H__
#define STATGEN_MI_GENOTYPE_DRAW_METHODS_H__

#include "statgen-mi/boost_random.h"

namespace MI {
  inline unsigned standard_genotype_draw(const double &p1, const double &p2) {
    double draw = MI::boost_random::next_real(0.0, 1.0);
    if (draw <= p1) return 1;
    if (draw <= p2) return 2;
    return 3;
  }
}

#endif //STATGEN_MI_GENOTYPE_DRAW_METHODS_H__
