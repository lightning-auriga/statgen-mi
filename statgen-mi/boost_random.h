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

#ifndef __MI_BOOST_RANDOM_H__
#define __MI_BOOST_RANDOM_H__
/*!
 * \file boost_random.h
 *
 * Interface for a better random number generator
 *
 */
#include "statgen-mi/config.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cfloat>
#include <ctime>
#include <cstdlib>
//#include <iostream>
#include <stdexcept>
namespace MI {
  /*!
   * \brief wrapper class for boost's rngs.
   *
   * This code is harvested from an integer rng implementation.
   * This wrapper class was originally designed to wrap
   * the variate generator in a method that is compatible with std::random_shuffle.
   * The class works just fine for the current purpose, however.  The call is
   *
   * boost_random::seed();
   * .....do other stuff.....
   * double my_random_number = boost_random::next_double();
   *
   * The length of the rand48 rng cycle is 2^48 - 1.  This is not
   * nearly as awesome as that of mt19937 (2^19937 - 1) but results
   * in a significant speed boost for the simulations; and while there
   * are intermediate distributions, they seem to not be present in all
   * versions of boost (ie taus88), which has caused some compiler issues.
   */
  class boost_random {
  public:
    /*!
     * \brief convenience typedef, per original recommendation on boost site
     */
    typedef boost::mt19937 base_generator_type;
    //typedef boost::rand48 base_generator_type;
    /*!
     * \brief convenience typedef, per original recommendation on boost site
     */
    typedef boost::uniform_int<> distribution_type;
    /*!
     * \brief convenience typedef, per original recommendation on boost site
     */
    typedef boost::variate_generator<base_generator_type&, distribution_type>
      generator_type;
    /*!
     * \brief static random number generator being wrapped by this class
     */
    static base_generator_type generator;
    /*!
     * \brief don't use this!  this is a static class really
     */
    boost_random() throw (std::domain_error) {
      throw std::domain_error("invalid random constructor");
    }
    /*!
     * \brief destructor: just for throughness
     */
    ~boost_random() throw () {
    }
    /*!
     * \brief seed based on system clock
     */
    static void seed() throw () {
      boost_random::generator.seed(static_cast<int32_t>(time(0)));
    }

    static double next_real(double min, double max) throw() {
      boost::variate_generator<base_generator_type &, boost::random::uniform_real_distribution<> > real_generator(boost_random::generator,
														  boost::random::uniform_real_distribution<>(min, max));
      double return_value = static_cast<double>(real_generator());
      return return_value;
    }

    /*!
     * \brief get the next random double on [0,1]
     * \return the next random double on [0,1]
     */
    static int next(std::ptrdiff_t arg) throw () {
      generator_type int_generator(boost_random::generator,
				   distribution_type(0, 
						     static_cast<unsigned>
						     (arg) - 1));
      int return_value =  static_cast<std::ptrdiff_t> (int_generator());
      return return_value;
    }
    static double next_normal(const double &mean, 
			      const double &variance) throw() {
      //std::cout << "running normal variate generator with mean " 
      //<< mean << " and variance " << variance << std::endl;
      boost::variate_generator<base_generator_type &,
	boost::random::normal_distribution<double> > gen
	(boost_random::generator, 
	 boost::random::normal_distribution<double>(mean,
						    sqrt(variance)));
      double return_value = static_cast<double>(gen());
      return return_value;
    }
    static unsigned next_geometric(const double &failp) throw() {
      boost::variate_generator<base_generator_type &,
	boost::random::geometric_distribution<unsigned, double> > gen
	(boost_random::generator,
	 boost::random::geometric_distribution<unsigned, double>(failp));
      unsigned return_value = static_cast<double>(gen());
      return return_value;
    }
  private:
  };
}
#endif // __MI_BOOST_RANDOM_H__
