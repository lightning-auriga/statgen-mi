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

#ifndef __MI_FILEINTERFACE_FILEINTERFACE_READER_FLAT_H__
#define __MI_FILEINTERFACE_FILEINTERFACE_READER_FLAT_H__

#include "fileinterface/fileinterface_reader_parent.h"

#include <string>
#include <fstream>
#include <stdexcept>

namespace MI {
  class fileinterface_reader_flat : public fileinterface_reader {
  public:
    fileinterface_reader_flat()
      : fileinterface_reader() {}
    ~fileinterface_reader_flat() throw() {}
    
    void open(const char *filename);
    void close() {_input.close(); clear();}
    void clear() {_input.clear();}
    bool is_open() const {return _input.is_open();}
    char get() {return _input.get();}
    bool getline(std::string &);
    bool eof() const {return _input.eof();}
    bool good() const {return _input.good();}
    bool bad() const {return _input.bad();}
    void read(char *, std::streamsize);
  private:
    std::ifstream _input;
  };
}

#endif //__MI_FILEINTERFACE_FILEINTERFACE_READER_FLAT_H__
