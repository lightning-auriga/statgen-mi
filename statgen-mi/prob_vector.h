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

/*!
 * \file prob_vector.h
 * \brief stored imputed data in compressed formate
 */
#ifndef STATGEN_MI_PROB_VECTOR_H__
#define STATGEN_MI_PROB_VECTOR_H__
#include <vector>
#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions/round.hpp>
namespace MI {
  /*!
   * \brief use stl vector to store data, but store it compressed to save lots of space
   */
  class prob_vector {
  public:
    /*!
     * \brief stl-like typedef for internal storage type
     */
    typedef unsigned value_type;
    /*!
     * \brief stl-like typedef for internal storage vector
     */
    typedef std::vector<value_type> vector;
    /*!
     * \brief stl-like typedef for internal storage size
     */
    typedef vector::size_type size_type;
    /*!
     * \brief constructor: initialize vector
     * @param init optional initial vector allocation size
     * @param n_bits the number of bits you want to use per value
     * \warning do NOT change n_bits unless you know what you're doing
     *
     * Setting the initial size can dramatically improve data insertion time
     */
    prob_vector(const std::vector<value_type>::size_type &init = 0,
		const size_type &n_bits = 10) throw(std::bad_alloc) 
      : _used_bound(0), _scaling(n_bits), _max_prec(1), _init(init) {
      reserve(init);
      _max_prec = (_max_prec << _scaling) - 1;
    }
    /*!
     * \brief copy constructor
     * @param other preexisting prob_vector
     */
    prob_vector(const prob_vector &other) throw(std::bad_alloc)
      : _data(other._data), _used_bound(other._used_bound), _scaling(other._scaling), _max_prec(other._max_prec) {}
    /*!
     * \brief destructor: all memory is managed
     */
    ~prob_vector() throw() {}




    /*!
     * \brief append a value to the end of the vector
     * @param value new datum to be added to vector
     * \warning any precision beyond [01].\d\d\d will be ignored
     */
    void push_back(const double &value) throw(std::bad_alloc);
    /*!
     * \brief set a value at a given index of the vector
     * @param index location in vector for insertion
     * @param value datum to be inserted
     */
    void set(const size_type &index, const double &value) throw(std::domain_error);
    /*!
     * \brief access a value at a given location in the vector
     * @param index location in the vector
     * \return stored value
     * \warning will throw std::domain_error if the index is invalid
     */
    double at(const size_type &index) const throw(std::domain_error);
    /*!
     * \brief determine how many values are stored
     * \return the number of stored values
     */
    inline size_type size() const throw() {return _used_bound;}


    /*!
     * \brief allocate but do not fill some data
     * @param amt amount of space to allocate
     *
     * this will not change data currently stored in the vector, but may require a copy
     */
    inline void reserve(const size_type &amt) throw(std::bad_alloc) {
      _data.reserve(static_cast<size_type>(ceil(amt * static_cast<double>(_scaling) / 8.0 / sizeof(value_type))));
    }

    inline void clear() {
      _used_bound = 0;
      _data.clear();
      reserve(_init);
    }
  private:
    /*!
     * \brief data vector
     */
    vector _data;
    /*!
     * \brief internal max_loc+1 of stored data
     */
    size_type _used_bound;
    /*!
     * \brief scaling factor for storage
     *
     * Basically, the number of bits you want
     */
    unsigned _scaling;
    /*!
     * \brief calculated maximum precision
     */
    unsigned _max_prec;
    /*!
     *
     */
    unsigned _init;
  };
}
#endif //STATGEN_MI_PROB_VECTOR_H__
