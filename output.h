/*
* Normaliz 2.2
* Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
* With contributions by Christof Soeger
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

#include "full_cone.h"
#include "lineare_transformation.h"
#include "sublattice_representation.h"


//---------------------------------------------------------------------------
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
  Full_Cone Result;
  Sublattice_Representation Basis_Change;
  bool BC_set;
//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//						Construction and destruction
//---------------------------------------------------------------------------

  Output();  //main constructor
  Output(const Output& Out);  //copy constructor
  ~Output();           		 //destructor

//---------------------------------------------------------------------------
//								Data acces
//---------------------------------------------------------------------------

  void read() const;                   // to be modified, just for tests
  void set_name(const string& n);             	//set name
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
  void set_result(const Full_Cone& C);         //sets Result
  void set_basis_change(const Sublattice_Representation& SR); // sets Basis_Change
  
  void write_matrix_ext(const Matrix& M) const; //writes M to file name.ext
  void write_matrix_ext_1(const Matrix& M) const; //writes M with a column of 1 added to file name.ext
  void write_matrix_esp(const Matrix& M) const; //writes M to file name.esp
  void write_matrix_typ(const Matrix& M) const; //writes M to file name.typ
  void write_matrix_egn(const Matrix& M) const; //writes M to file name.egn
  void write_matrix_gen(const Matrix& M) const; //writes M to file name.gen
  void write_matrix_sup(const Matrix& M) const; //writes M to file name.supwrite_out
  void write_matrix_tri(const Matrix& M) const; //writes M to file name.tri
  void write_matrix_ht1(const Matrix& M) const; //writes M to file name.tri

//---------------------------------------------------------------------------
//                         Output Algorithms
//---------------------------------------------------------------------------

  void cone()const;
  void polytop()const;
  void rees(const bool primary)const;
  void dual()const;

//---------------------------------------------------------------------------
//							Error msg
//---------------------------------------------------------------------------

  void error(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

