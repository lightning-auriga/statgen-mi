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

#ifndef __MI_TRANSPOSE_MATRIX_H__
#define __MI_TRANSPOSE_MATRIX_H__

#include <vector>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "statgen-mi/prob_vector.h"
#include "statgen-mi/utilities.h"
#include "statgen-mi/annotations.h"
#include "fileinterface/fileinterface.h"

namespace MI {
  class matrix_transposer {
  public:
#ifdef FILEINTERFACE_HAVE_LIBZ
    typedef MI::fileinterface_reader_gzip fileinterface_reader_type;
#else
    typedef MI::fileinterface_reader_flat fileinterface_reader_type;
#endif //FILEINTERFACE_HAVE_LIBZ
    typedef typename std::vector<std::pair<std::string, boost::shared_ptr<fileinterface_reader_type> > > subunit_vector;
    matrix_transposer()
      : _filename(""),
      _input(0),
      _current_column(0) {
      throw std::domain_error("MI::matrix_transposer: must be initialized with filename");
    }
    matrix_transposer(const std::string &filename)
      : _filename(filename),
      _input(0),
      _current_column(0) {}

    ~matrix_transposer() throw() {clear_internals();}

    bool transpose_a_matrix(void (*input_converter)(const std::string &, MI::prob_vector &, MI::annotations &),
			    //void (*output_converter)(const MI::prob_vector &, std::string &),
			    prob_vector &result,
			    unsigned ram_limit = 4,
			    unsigned disk_limit = 20);

    bool get_a_line(void (*input_converter)(const std::string &, MI::prob_vector &, MI::annotations &),
		    //void (*output_converter)(const MI::prob_vector &, std::string &),
		    prob_vector &result,
		    annotations &annot,
		    bool input_is_snp_major,
		    bool output_is_snp_major);

    void clear_internals();

    const std::string &get_filename() const {return _filename;}
  private:
    void write_current_matrix();
    bool get_a_line_from_single_matrix(prob_vector &result);
    bool get_a_line_from_multiple_files(prob_vector &result);
    
    std::string _filename;
    MI::fileinterface_reader *_input;
    subunit_vector _matrix_subunits;
    std::vector<MI::prob_vector> _single_matrix;
    unsigned _current_column;

    std::vector<MI::annotations> _per_line_annotations;
  };
}

#endif //__MI_TRANSPOSE_MATRIX_H__
