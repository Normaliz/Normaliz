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

//---------------------------------------------------------------------------
#ifndef MAP_OPERATIONS_H
#define MAP_OPERATIONS_H


//---------------------------------------------------------------------------
                  
#include <map>
#include <ostream>

namespace libnormaliz {
using std::ostream;

template<typename key, typename T>
ostream& operator<< (ostream& out, const map<key, T> M) {
    typename map<key, T>::const_iterator it;
    for (it = M.begin(); it != M.end(); ++it) {
        out << it->first << ": " << it-> second << "  ";
    }
    out << std::endl;
    return out;
}

}  //end namespace

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
