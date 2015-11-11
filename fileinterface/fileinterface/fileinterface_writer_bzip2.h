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

#ifndef __MI_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__
#define __MI_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__

#include "fileinterface/config.h"
#ifdef FILEINTERFACE_HAVE_LIBBZ2

#include "fileinterface/fileinterface_writer_parent.h"

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <bzlib.h>

namespace MI {
  class fileinterface_writer_bzip2 : public fileinterface_writer {
  public:
    fileinterface_writer_bzip2()
      : fileinterface_writer(),
      _raw_output(0), _bz_output(0) {}
    ~fileinterface_writer_bzip2() throw() {close();}

    void open(const char *filename);
    void close();
    void clear();
    bool is_open() const;
    void put(char c);
    void writeline(const std::string &);
    void write(char *, std::streamsize);
    bool eof() const {return false;}
    bool good() const {return _good;}
    bool fail() const {return _fail;}
    bool bad() const {return _bad;}
  private:
    FILE *_raw_output;
    BZFILE *_bz_output;
  };
}


#endif //FILEINTERFACE_HAVE_LIBBZ2

#endif //__MI_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__
