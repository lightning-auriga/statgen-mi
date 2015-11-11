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

#ifndef __MI_FILEINTERFACE_FILEINTERFACE_READER_PARENT_H__
#define __MI_FILEINTERFACE_FILEINTERFACE_READER_PARENT_H__

#include <string>
#include <ios>

#include "fileinterface/helper.h"

namespace MI {
  class fileinterface_reader {
  public:
  fileinterface_reader() 
    : _good(true), _bad(false), _fail(false) {}
    ~fileinterface_reader() throw() {}

    void open(const std::string &filename) {open(filename.c_str());}
    virtual void open(const char *filename) = 0;
    virtual void close() = 0;
    virtual void clear() = 0;
    virtual bool is_open() const = 0;
    virtual char get() = 0;
    virtual bool getline(std::string &) = 0;
    virtual bool eof() const = 0;
    virtual bool good() const = 0;
    virtual bool bad() const = 0;
    virtual void read(char *, std::streamsize) = 0;
  protected:
    bool _good;
    bool _bad;
    bool _fail;
  };
}

#endif //__MI_FILEINTERFACE_FILEINTERFACE_READER_PARENT_H__
