/*
 * Normaliz 2.5
 * Copyright (C) 2007-2010  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
#ifndef OUTPUT_H
#define OUTPUT_H
//---------------------------------------------------------------------------

#include "libnormaliz/integer.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/lineare_transformation.h"
#include "libnormaliz/sublattice_representation.h"
#include "libnormaliz/cone.h"

using namespace std;
using namespace libnormaliz;

//---------------------------------------------------------------------------

template<typename Integer>
class Output {
  string name;
  bool out;
  bool inv;
  bool ext;
  bool esp;
  bool typ;
  bool egn;
  bool gen;
  bool sup;
  bool tri;
  bool ht1;
  Cone<Integer>* Result;

//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//						Construction and destruction
//---------------------------------------------------------------------------

  Output();  //main constructor
  Output(const Output<Integer>& Out);  //copy constructor
  ~Output();           		 //destructor

//---------------------------------------------------------------------------
//								Data acces
//---------------------------------------------------------------------------

  void read() const;                   // to be modified, just for tests
  void set_name(const string& n);             	//set name
  void setCone(Cone<Integer> & C);
  void set_write_out(const bool& flag);             //sets the write .out flag
  void set_write_inv(const bool& flag);             //sets the write .inv flag
  void set_write_ext(const bool& flag);             //sets the write .ext flag
  void set_write_esp(const bool& flag);             //sets the write .esp flag
  void set_write_typ(const bool& flag);             //sets the write .typ flag
  void set_write_egn(const bool& flag);             //sets the write .egn flag
  void set_write_gen(const bool& flag);             //sets the write .gen flag
  void set_write_sup(const bool& flag);             //sets the write .sup flag
  void set_write_tri(const bool& flag);             //sets the write .tri flag
  void set_write_ht1(const bool& flag);             //sets the write .ht1 flag
  void set_write_extra_files();         	    //sets some flags to true
  void set_write_all_files();          		    //sets all flags to true
  
  void write_matrix_ext(const Matrix<Integer>& M) const; //writes M to file name.ext
  void write_matrix_esp(const Matrix<Integer>& M) const; //writes M to file name.esp
  void write_matrix_typ(const Matrix<Integer>& M) const; //writes M to file name.typ
  void write_matrix_egn(const Matrix<Integer>& M) const; //writes M to file name.egn
  void write_matrix_gen(const Matrix<Integer>& M) const; //writes M to file name.gen
  void write_matrix_sup(const Matrix<Integer>& M) const; //writes M to file name.sup
  void write_matrix_tri(const Matrix<Integer>& M) const; //writes M to file name.tri
  void write_matrix_ht1(const Matrix<Integer>& M) const; //writes M to file name.ht1

  void write_inv_file() const;


//---------------------------------------------------------------------------
//                         Output Algorithms
//---------------------------------------------------------------------------

  void cone()const;
  void polytop()const;
  void rees(const bool primary)const;

};
//class end *****************************************************************

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

