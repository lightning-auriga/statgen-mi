/*
  Copyright 2022 Lightning Auriga

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

#ifndef STATGEN_MI_DESCRIPTIONS_H__
#define STATGEN_MI_DESCRIPTIONS_H__

#include <string>
#include <stdexcept>

namespace MI {
  class program_description {
  public:
  program_description()
    : _snp_major(false),
      _program_exec(""),
      _program_handle("") {}
    ~program_description() throw() {}
    program_description(const program_description &obj)
      : _snp_major(obj._snp_major),
      _program_exec(obj._program_exec) {}

    inline bool is_snp_major() const {return _snp_major;}
    inline bool is_indiv_major() const {return !is_snp_major();}
    inline const std::string &program_executable() const {return _program_exec;}
    inline const std::string &program_handle() const {return _program_handle;}

    inline void set_snp_major(bool b) {_snp_major = b;}
    inline void set_program_executable(const std::string &s) {_program_exec = s;}
    inline void set_program_handle(const std::string &s) {_program_handle = s;}
  private:
    bool _snp_major;
    std::string _program_exec;
    std::string _program_handle;
  };

  class imputed_dataset_description {
  public:
    imputed_dataset_description()
      : _snp_major(false),
      _program_handle("") {}
    ~imputed_dataset_description() throw() {}
    imputed_dataset_description(const imputed_dataset_description &obj)
      : _snp_major(obj._snp_major),
      _program_handle(obj._program_handle) {}

    inline bool is_snp_major() const {return _snp_major;}
    inline bool is_indiv_major() const {return !is_snp_major();}
    inline const std::string &program_handle() const {return _program_handle;}

    inline void set_snp_major(bool b) {_snp_major = b;}
    inline void set_program_handle(const std::string &s) {_program_handle = s;}
  private:
    bool _snp_major;
    std::string _program_handle;
  };
}

#endif //STATGEN_MI_DESCRIPTIONS_H__
