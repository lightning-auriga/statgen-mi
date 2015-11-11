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

#include "statgen-mi/submission_formatter.h"

std::string MI::submission_formatter::get_submission_command(const std::string &queue_type,
							     const std::string &command,
							     const std::string &output_filename,
							     const std::string &error_filename,
							     const std::string &mem_limit,
							     const std::string &time_limit) const {
  unsigned mem_in_megabytes = reconcile_memory(mem_limit);
  unsigned time_in_hours = reconcile_time(time_limit);
  if (cicompare(queue_type, "bsub")) {
    return get_submission_command_bsub(command,
				       output_filename,
				       error_filename,
				       mem_in_megabytes,
				       time_in_hours);
  } else if (cicompare(queue_type, "qsub")) {
    return get_submission_command_qsub(command,
				       output_filename,
				       error_filename,
				       mem_in_megabytes,
				       time_in_hours);
  } else if (cicompare(queue_type, "none")) {
    return get_submission_command_none(command,
				       output_filename,
				       error_filename);
  } else {
    throw std::domain_error("submission_formatter: unrecognized queue type: \"" + queue_type + "\"");
  }
}

std::string MI::submission_formatter::get_submission_command_none(const std::string &command,
								  const std::string &output_filename,
								  const std::string &error_filename) const {
  return command + " && echo \"Successfully completed\" > " + output_filename;
}

std::string MI::submission_formatter::get_submission_command_bsub(const std::string &command,
								 const std::string &output_filename,
								 const std::string &error_filename,
								 unsigned           mem_limit,
								 unsigned           time_limit) const {
  std::string prefix = "bsub -q " + get_best_queue("bsub", mem_limit, time_limit)
    + " -o " + output_filename
    + " -e " + error_filename
    + " -R \"rusage[mem=" + to_string<unsigned>(mem_limit)
    + "M]\" -M " + to_string<unsigned>(mem_limit) + "M ";
  return prefix + command;
}

std::string MI::submission_formatter::get_submission_command_qsub(const std::string &command,
								  const std::string &output_filename,
								  const std::string &error_filename,
								  unsigned           mem_limit,
								  unsigned           time_limit) const {
  std::string suffix = "qsub -l mem=" + to_string<unsigned>(mem_limit)
    + "M,time=" + to_string<unsigned>(time_limit) + ":: -cwd ";
  return "echo `" + command + " && echo \\\"Successfully completed\\\" > " + output_filename + "` | " + suffix;
}








unsigned MI::submission_formatter::reconcile_memory(const std::string &mem) const {
  //want memory in MB
  /* accept formats:
     ### or ####: a number in MB
     # or ##: a number in GB
     (#+)[GM]: a number in the specified unit
   */
  unsigned res = 0;
  std::string running_mem = mem;
  if (running_mem.find_first_not_of("0123456789") == std::string::npos) {
    res = from_string<unsigned>(running_mem);
    if (res < 100) res = res * 1000;
  } else if (running_mem.find("M") != std::string::npos) {
    running_mem = running_mem.substr(0, running_mem.find("M"));
    res = from_string<unsigned>(running_mem);
  } else if (running_mem.find("G") != std::string::npos) {
    running_mem = running_mem.substr(0, running_mem.find("G"));
    res = from_string<unsigned>(running_mem) * 1000;
  } else {
    throw std::domain_error("submission_formatter: unrecognized memory format: \"" + mem + "\"");
  }
  return res;
}

unsigned MI::submission_formatter::reconcile_time(const std::string &time) const {
  //want time in hours
  /* accept formats:
     (#+): a count of hours
     (#+):(#*):(#*) | (#*):(#+):(#*) | (#*):(#*):(#+): HH:MM:SS basically
   */
  unsigned res = 0;
  if (time.find_first_not_of("0123456789") == std::string::npos) {
    res = from_string<unsigned>(time);
  } else if (time.find(":") != std::string::npos) {
    std::string shours = "", sminutes = "", sseconds = "", adj_time = "";
    unsigned hours = 0, minutes = 0, seconds = 0;
    shours = time.substr(0, time.find(":"));
    adj_time = time.substr(time.find(":") + 1);
    if (adj_time.find(":") != std::string::npos) {
      sminutes = adj_time.substr(0, adj_time.find(":"));
      sseconds = adj_time.substr(adj_time.find(":") + 1);

      if (shours.empty()) {
	hours = 0;
      } else {
	hours = from_string<unsigned>(shours);
      }
      if (sminutes.empty()) {
	minutes = 0;
      } else {
	minutes = from_string<unsigned>(sminutes);
      }
      if (sseconds.empty()) {
	seconds = 0;
      } else {
	seconds = from_string<unsigned>(sseconds);
      }
      res = hours + minutes / 60 + seconds / 3600;
    } else {
      throw std::domain_error("submission_formatter: partial time format detected, please correct: \"" + time + "\"");
    }
    return res;
  } else {
    throw std::domain_error("submission_formatter: unrecognized time format: \"" + time + "\"");
  }
}


void MI::submission_formatter::initialize_queue_maps() {
  std::string queue_config_filename = MI::parameters::get_parameter("mi-submission-queue-config");
  MI::fileinterface_reader *input = 0;
  std::string queue_name = "", queue_type = "", line = "";
  unsigned queue_max_mem = 0, queue_max_time = 0;
  try {
    input = MI::reconcile_reader(queue_config_filename);
    while (input->getline(line)) {
      if (line.find("#") != std::string::npos) line = line.substr(0, line.find("#"));
      if (line.empty()) continue;
      std::istringstream strm1(line);
      if (!(strm1 >> queue_name >> queue_type >> queue_max_mem >> queue_max_time))
	throw std::domain_error("submission_formatter: could not parse config line \"" + line + "\"");
      if (queue_memory_mapper.find(queue_type) == queue_memory_mapper.end()) {
	queue_memory_mapper[queue_type] = std::multimap<unsigned, std::string>();
      }
      if (queue_time_mapper.find(queue_type) == queue_time_mapper.end()) {
	queue_time_mapper[queue_type] = std::multimap<unsigned, std::string>();
      }
      queue_memory_mapper[queue_type].insert(std::make_pair(queue_max_mem, queue_name));
      queue_time_mapper[queue_type].insert(std::make_pair(queue_max_time, queue_name));
    }
    input->close();
    delete input;
  } catch (...) {
    if (input) delete input;
    //throw;
    std::cerr << "warning: no queue map provided; if you're using a submission mode that doesn't require this, you're fine" << std::endl;
  }
}

std::string MI::submission_formatter::get_best_queue(const std::string &queue_type,
						     unsigned           mem_limit,
						     unsigned           time_limit) const {
  std::map<std::string, std::multimap<unsigned, std::string> >::const_iterator mem_queue, time_queue;
  std::multimap<unsigned, std::string>::const_iterator mem_iter, time_iter;
  std::map<std::string, bool> permissible_queues;
  if ((mem_queue = queue_memory_mapper.find(queue_type)) == queue_memory_mapper.end() ||
      (time_queue = queue_time_mapper.find(queue_type)) == queue_time_mapper.end())
    throw std::domain_error("submission_formatter::get_best_queue: requested queue type \""
			    + queue_type + "\" has no stored information");
  mem_iter = mem_queue->second.lower_bound(mem_limit);
  for ( ; mem_iter != mem_queue->second.end(); ++mem_iter) {
    permissible_queues[mem_iter->second] = true;
  }
  time_iter = time_queue->second.lower_bound(time_limit);
  for ( ; time_iter != time_queue->second.end(); ++time_iter) {
    if (permissible_queues.find(time_iter->second) != permissible_queues.end())
      return time_iter->second;
  }
  throw std::domain_error("submission_formatter: unable to find permissible queue given provided "
			  "memory and time requirements");
}
