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

#ifndef __MI_FILEINTERFACE_FILEINTERFACE_WRITER_FLAT_H__
#define __MI_FILEINTERFACE_FILEINTERFACE_WRITER_FLAT_H__

#include <fstream>
#include <string>
#include <stdexcept>
#include "fileinterface/fileinterface_writer_parent.h"

namespace MI {
  class fileinterface_writer_flat : public fileinterface_writer {
  public:
    fileinterface_writer_flat()
      : fileinterface_writer() {}
    ~fileinterface_writer_flat() throw() {close();}

    void open(const char *filename);
    void close() {_output.close(); clear();}
    void clear() {_output.clear();}
    bool is_open() const {return _output.is_open();}
    void put(char c) {_output.put(c);}
    void writeline(const std::string &);
    void write(char *, std::streamsize);
    bool eof() const {return _output.eof();}
    bool good() const {return _output.good();}
    bool fail() const {return _output.fail();}
    bool bad() const {return _output.bad();}
  private:
    std::ofstream _output;
  };
}

#endif //__MI_FILEINTERFACE_FILEINTERFACE_WRITER_FLAT_H__
