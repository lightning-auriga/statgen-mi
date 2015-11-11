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

#include "fileinterface/fileinterface_writer_flat.h"

void MI::fileinterface_writer_flat::open(const char *filename) {
  _output.open(filename);
  if (!_output.is_open())
    throw std::domain_error("MI::fileinterface_writer_flat::open: cannot "
			    "write file \"" + std::string(filename) + "\"");
}

void MI::fileinterface_writer_flat::writeline(const std::string &line) {
  _output << line << get_newline();
}

void MI::fileinterface_writer_flat::write(char *buf, std::streamsize n) {
  _output.write(buf, n);
}
