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

#ifndef __MI_HELPER_H__
#define __MI_HELPER_H__

#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <boost/algorithm/string/find.hpp>


namespace MI {
  template <class value_type>
    std::string to_string(const value_type &obj) {
    std::ostringstream o;
    if (!(o << obj))
      throw std::domain_error("to_string: cannot convert object to string");
    return o.str();
  }

  template <class value_type>
    value_type from_string(const std::string &str) {
    std::istringstream strm1(str);
    value_type res;
    if (!(strm1 >> res))
      throw std::domain_error("from_string: cannot convert string to object: \""
			      + str + "\"");
    return res;
  }

  inline std::string::const_iterator find_nth_occurrence(const std::string &str,
							 const std::string &substr,
							 unsigned           which) {
    boost::iterator_range<std::string::const_iterator> iter = boost::algorithm::find_nth(str, substr, which);
    return iter.begin();
  }

  inline std::string substring_after_nth_occurrence(const std::string &str,
						    const std::string &substr,
						    unsigned           which) {
    std::string::const_iterator iter = find_nth_occurrence(str, substr, which);
    return std::string(iter + 1, str.end());
  }

  inline std::string substring_before_nth_occurrence(const std::string &str,
						     const std::string &substr,
						     unsigned           which) {
    std::string::const_iterator iter = find_nth_occurrence(str, substr, which);
    return std::string(str.begin(), iter);
  }

  inline unsigned max(unsigned u1, unsigned u2) throw() {
    return u1 > u2 ? u1 : u2;
  }

  inline double max(const double &d1, const double &d2) throw() {
    return d1 > d2 ? d1 : d2;
  }

  inline unsigned convert_to_tabs(unsigned filled_chars, unsigned max_tabs, unsigned space_per_tab) {
    unsigned filled_tabs = filled_chars / space_per_tab;
    return max_tabs - filled_tabs;
  }

  inline bool cicompare(const std::string &str1, const std::string &str2) {
    if (str1.size() == str2.size()) {
      for (unsigned i = 0; i < str1.size(); ++i) {
	if (tolower(str1.at(i)) != tolower(str2.at(i))) return false;
      }
      return true;
    }
    return false;
  }

  inline bool pair_double_first_sort(const std::pair<double, unsigned> &p1,
				     const std::pair<double, unsigned> &p2) {
    return p1.first > p2.first;
  }

  inline double logit(const double &d) {
    return log(d / (1.0-d));
  }

  inline double unlogit(const double &d) {
    return exp(d) / (1.0 + exp(d));
  }

  inline std::pair<double, double> from_y_z(const double &y, const double &z) {
    std::pair<double, double> res;
    res.first = z*y/(1.0-y+y*z);
    res.second = z*(1.0-y)/(1.0-y+y*z); //, 1.0-z/(1.0-y+y*z);
    return res;
  }

  inline std::pair<double, double> to_y_z(const double &p1, const double &p2) {
    std::pair<double, double> res;
    res.first = p1/(p1 + p2);
    res.second = p2/(1.0 - p1);
    return res;
  }

  inline double process_unit_interval(const double &d, const double &unit_factor) {
    return d*(1.0 - unit_factor) + unit_factor / 2.0;
  }

  inline std::pair<double, double> deprocess_pvalues(const double &p1, const double &p2, const double &unit_factor) {
    std::pair<double, double> res;
    res.first = (p1 - unit_factor/3.0) / (1.0 - unit_factor);
    res.second = (p2 - unit_factor/3.0) / (1.0 - unit_factor);
    return res;
  }

  inline std::pair<double, double> process_pvalues(const double &p1, const double &p2, const double &unit_factor) {
    std::pair<double, double> res;
    res.first = p1 * (1.0 - unit_factor) + unit_factor / 3.0;
    res.second = p2 * (1.0 - unit_factor) + unit_factor / 3.0;
    return res;
  }

  inline void string_to_vector(const std::string &line, std::vector<std::string> &vec, char sep = '\t') {
    vec.clear();
    std::string substring = "";
    std::string::size_type loc = 0, cur = 0;
    if (line.find(sep) == std::string::npos) {
      vec.push_back(line);
    } else {
      while ((loc = line.find(sep)) != std::string::npos) {
	substring = line.substr(cur, loc);
	vec.push_back(substring);
	cur = loc;
      }
    }
  }

  inline bool successfully_completed(const std::string &filename) {
    std::ifstream input;
    std::string line = "";
    input.open(filename.c_str());
    if (!input.is_open()) return false;
    while (input.peek() != EOF) {
      getline(input, line);
      if (line.find("Successfully completed") != std::string::npos) {
	input.close();
	return true;
      }
    }
    input.close();
    return false;
  }

  inline void safely_remove(const std::string &filename,
			    bool error_on_fail = true) {
    std::istringstream strm1(filename);
    std::string indiv_filename = "";
    while (strm1 >> indiv_filename) {
      if (remove(indiv_filename.c_str()) && error_on_fail)
	throw std::domain_error("unable to remove file \"" + filename + "\"");
    }
  }
}
#endif //__MI_HELPER_H__
