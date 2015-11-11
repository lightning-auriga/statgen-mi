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

#ifndef __MI_SUBMISSION_FORMATTER_H__
#define __MI_SUBMISSION_FORMATTER_H__

#include <string>
#include <stdexcept>
#include <map>
#include "statgen-mi/utilities.h"
#include "statgen-mi/cargs.h"
#include "fileinterface/fileinterface.h"

namespace MI {

class submission_formatter {
 public:
  submission_formatter() {initialize_queue_maps();}
  ~submission_formatter() throw() {}

  //within this method, add custom queue handling from queue_type
  std::string get_submission_command(const std::string &queue_type,
				     const std::string &command,
				     const std::string &output_filename,
				     const std::string &error_filename,
				     const std::string &mem_limit,
				     const std::string &time_limit) const;

  ///////////// and then add your own queue handler down here
  std::string get_submission_command_bsub(const std::string &command,
					  const std::string &output_filename,
					  const std::string &error_filename,
					  unsigned           mem_limit,
					  unsigned           time_limit) const;
  std::string get_submission_command_qsub(const std::string &command,
					  const std::string &output_filename,
					  const std::string &error_filename,
					  unsigned           mem_limit,
					  unsigned           time_limit) const;
  std::string get_submission_command_none(const std::string &command,
					  const std::string &output_filename,
					  const std::string &error_filename) const;
  ///////////////////////////////////////////////////////////
private:
  std::map<std::string, std::multimap<unsigned, std::string> > queue_memory_mapper;
  std::map<std::string, std::multimap<unsigned, std::string> > queue_time_mapper;

  void initialize_queue_maps();
  unsigned reconcile_memory(const std::string &str) const;
  unsigned reconcile_time(const std::string &str) const;
  std::string get_best_queue(const std::string &queue_type,
			     unsigned           mem_limit,
			     unsigned           time_limit) const;
};

}
#endif //__MI_SUBMISSION_FORMATTER_H__
