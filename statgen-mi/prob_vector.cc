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

#include "statgen-mi/prob_vector.h"

void statgen_mi::prob_vector::push_back(const double &value) throw(std::bad_alloc) {
  //assuming 10bits fits in a single size_type
  if (size() >= floor(8.0 / static_cast<double>(_scaling) * _data.size() * sizeof(value_type)))
    _data.push_back(0);
  ++_used_bound;
  set(size()-1, value);
}

void statgen_mi::prob_vector::set(const prob_vector::size_type &index, const double &value) throw(std::domain_error) {
  if (index >= size()) throw std::domain_error("prob_vector set: index out of range");
  value_type converted = static_cast<value_type>(boost::math::lround<double>(1000.0 * value));
  if (converted < 0 || converted > (_max_prec-1)) converted = _max_prec - 1;
  size_type x = (_scaling * index) % (sizeof(value_type) << 3);
  size_type y = (_scaling * index) / (sizeof(value_type) << 3);
  size_type z = (_scaling + x >= (sizeof(value_type) << 3) ? (_scaling + x) - (sizeof(value_type) << 3) : 0);
  value_type converted_shift = ((sizeof(value_type) << 3) - x < _scaling ?
				(converted >> (_scaling - (sizeof(value_type) << 3) + x)) :
				(converted << ((sizeof(value_type) << 3) - x - _scaling)));
  value_type mask = (1 << _scaling) - 1;
  mask = ((sizeof(value_type) << 3) - x < _scaling ?
	  (mask >> (_scaling - (sizeof(value_type) << 3) + x)) :
	  (mask << ((sizeof(value_type) << 3) - x - _scaling)));
  _data.at(y) = _data.at(y) & ~mask;
  _data.at(y) = _data.at(y) | converted_shift;
  if (z > 0) {
    _data.at(y + 1) = _data.at(y + 1) | (converted << ((sizeof(value_type) << 3) - z));
  }
}

double statgen_mi::prob_vector::at(const prob_vector::size_type &index) const throw(std::domain_error) {
  if (index >= size()) throw std::domain_error("prob_vector at: index out of range");
  value_type converted = (1 << _scaling) - 1;
  size_type x = (_scaling * index) % (sizeof(value_type) << 3);
  size_type y = (_scaling * index) / (sizeof(value_type) << 3);
  size_type z = (_scaling + x >= (sizeof(value_type) << 3) ? (_scaling + x) - (sizeof(value_type) << 3) : 0);
  converted = converted & ((sizeof(value_type) << 3) - x < _scaling ?
			   (_data.at(y) << (_scaling - (sizeof(value_type) << 3) + x)) :
			   (_data.at(y) >> ((sizeof(value_type) << 3) - x - _scaling)));
  if (z > 0) {
    converted = converted | (_data.at(y + 1) >> ((sizeof(value_type) << 3) - z));
  }
  return converted > _max_prec - 1 ? -9.0 : static_cast<double>(converted)/1000.0;
}
