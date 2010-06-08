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
#ifndef MODE_H
#define MODE_H
//---------------------------------------------------------------------------
#include "output.h"
#include "matrix.h"

//---------------------------------------------------------------------------
//                              Mode selection
//---------------------------------------------------------------------------

void make_main_computation(const int& mode, string& run_mode_type, const Matrix& Input, Output& Out);
void run_mode_0( string& run_mode_type,const Matrix& Input, Output& Out);
void run_mode_1( string& run_mode_type,const Matrix& Input, Output& Out);
void run_mode_2( string& run_mode_type,const Matrix& Input, Output& Out);
void run_mode_3( string& run_mode_type,const Matrix& Input, Output& Out);
void run_mode_10( string& computation_type,const Matrix& Binomials, Output& Out);

void run_mode_456(string& computation_type, const Matrix& Congruences, Matrix Equations, Matrix Inequalities, Output& Out);
void run_mode_4( string& run_mode_type,const Matrix& Input, const int& nr_equations, Output& Out);
void run_mode_5( string& run_mode_type,const Matrix& Input, Output& Out);
void run_mode_equ_inequ( string& computation_type,const Matrix& Equations, const Matrix& Inequalities, Output& Out);

//---------------------------------------------------------------------------
//                          Run_mode_type selection
//---------------------------------------------------------------------------

Full_Cone make_computations(const string& run_mode_type, const Matrix& Full_Cone_Generators);

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
