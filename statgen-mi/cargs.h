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

/*!
 * \file cargs.h
 * \brief contains argument parsing class declarations
 */

#ifndef __MI_CARGS_H__
#define __MI_CARGS_H__

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>
#include <stdexcept>

#include "statgen-mi/config.h"
#include "statgen-mi/parameters.h"
#include "statgen-mi/utilities.h"
//#include "statgen-mi/program_configuration.h"

namespace MI {
  class cargs {
  public:
    /*!
     * \brief default constructor; don't really use this...
     */
    cargs() throw(std::domain_error) {
      throw std::domain_error("Programmer: do not use default constructor of cargs");
    }
    /*!
     * \brief useful constructor;
     * @param argc from main, number of command line arguments, including program call
     * @param argv from main, pointer to array of command line arguments
     */
    cargs(int argc, char **argv) throw(std::domain_error) :
    _argc(argc), _argv(argv) {
      set_defaults();
    }
    /*!
     * \brief boring destructor
     */
    ~cargs() throw() {
    }
    /*!
     * \brief call after setting argc and argv.  Parses command line arguments
     */
    void parse_args(bool warn_unused) throw(std::domain_error);
    void print_ignored() const throw();
    void set_flag(const std::string &parameter_long, const std::string &parameter_short, const std::string &description, bool default_value);
    void set_parameter(const std::string &parameter_long, const std::string &parameter_short, const std::string &description, const std::string &default_value = "");

  protected:
    /*!
     * \brief access argc
     * @return argc
     */
    int num_args() const throw() {
      return _argc;
    }
    /*!
     * \brief access argv
     * @return argv
     */
    char **arg_array() const throw() {
      return _argv;
    }
  private:
    /*!
     * \brief initialize some private data members
     */
    void set_defaults() throw(std::domain_error);
    /*!
     * \brief search the command line parameters for a given flag
     * @param key search term
     * @return whether or not the term was found
     */
    bool find(const std::string &key) throw(std::domain_error);
    /*!
     * \brief once you confirm there is a key, get a value_type right after the key
     * @param key search term (already used in a find call)
     * @return formatted value_type following the key in the command line argument array
     */
    template<class value_type>
      value_type get_formatted_value(const std::string &key) throw(std::domain_error) {
      for (int i = 1; i < num_args(); ++i) {
	if (!std::string(arg_array()[i]).compare(key)) {
	  if (!_found.at(i)) {
	    throw std::domain_error("Programmer: must call find on key \"" + key
				    + "\" before getting value");
	  }
	  if (i == num_args() - 1) {
	    throw std::domain_error("User: command line flag \"" + key
				    + "\" expects an argument but has none");
	  }
	  _found.at(i + 1) = true;
	  return from_string<value_type> (std::string(arg_array()[i + 1]));
	}
      }
      throw std::domain_error("Programmer: consistency std::domain_error: cannot get value for key \"" + key + "\"");
    }
    int _argc;
    char **_argv;
    std::vector<bool> _found;
    //program_configuration _configuration;
  };
}
#endif /* __MI_CARGS_H__ */
