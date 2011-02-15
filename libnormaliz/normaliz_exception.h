/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef NORMALIZ_EXEPTION_H_
#define NORMALIZ_EXEPTION_H_

#include <exception>
#include "libnormaliz.h"

namespace libnormaliz {

class NormalizException: public virtual std::exception {
	virtual const char* what() const throw() {
		return "NormalizException happened!";
	}
};

class ArithmeticException: public virtual NormalizException {
	virtual const char* what() const throw() {
		return "Arithmetic Overflow detected, try a bigger integer type!";
	}
};



} /* end namespace */

#endif /* LIBNORMALIZ_H_ */
