/*
  Copyright 2022 Lightning Auriga
 
  This file is part of statgen-mi
 
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

#ifndef STATGEN_MI_ANNOTATIONS_H__
#define STATGEN_MI_ANNOTATIONS_H__

#include <map>
#include <string>
#include <vector>

namespace statgen_mi {
  class annotations {
  public:
    annotations() {}
    annotations(const annotations &obj)
      : _data(obj._data) {}
    ~annotations() throw() {}

    void set(const std::string &key, const std::string &val) {
      std::vector<std::string> vec;
      vec.push_back(val);
      _data[key] = vec;
    }
    void add(const std::string &key, const std::string &val) {
      annotation_struct::iterator iter = _data.find(key);
      if (iter == _data.end()) {
	set(key, val);
      } else {
	iter->second.push_back(val);
      }
    }
    std::vector<std::string> get(const std::string &key) const {
      annotation_struct::const_iterator iter;
      if ((iter = _data.find(key)) != _data.end()) return iter->second;
      return std::vector<std::string>();
    }
    std::string get(const std::string &key, unsigned index) const {
      annotation_struct::const_iterator iter = _data.find(key);
      if (iter != _data.end()) {
	if (index >= iter->second.size()) return "";
	return iter->second.at(index);
      } else return "";
    }
    unsigned get_count_at(const std::string &key) const {
      annotation_struct::const_iterator iter = _data.find(key);
      if (iter == _data.end()) return 0;
      return iter->second.size();
    }
    void clear() {_data.clear();}
  private:
    typedef std::map<std::string, std::vector<std::string> > annotation_struct;
    annotation_struct _data;
  };
}

#endif //STATGEN_MI_ANNOTATIONS_H__
