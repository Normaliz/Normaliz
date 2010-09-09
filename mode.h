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
#ifndef MODE_H
#define MODE_H
//---------------------------------------------------------------------------
#include "libnormaliz.h"
#include "output.h"
#include "matrix.h"

//---------------------------------------------------------------------------
//                              Mode selection
//---------------------------------------------------------------------------

template<typename Integer>
void make_main_computation(const int& mode, string& run_mode_type, const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_0( string& run_mode_type,const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_1( string& run_mode_type,const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_2( string& run_mode_type,const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_3( string& run_mode_type,const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_10( string& computation_type,const Matrix<Integer>& Binomials, Output<Integer>& Out);

template<typename Integer>
void run_mode_456(string& computation_type, const Matrix<Integer>& Congruences, Matrix<Integer> Equations, Matrix<Integer> Inequalities, Output<Integer>& Out);
template<typename Integer>
void run_mode_4( string& run_mode_type,const Matrix<Integer>& Input, const int& nr_equations, Output<Integer>& Out);
template<typename Integer>
void run_mode_5( string& run_mode_type,const Matrix<Integer>& Input, Output<Integer>& Out);
template<typename Integer>
void run_mode_equ_inequ( string& computation_type,const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities, Output<Integer>& Out);

//---------------------------------------------------------------------------
//                          Run_mode_type selection
//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer> make_computations(const string& run_mode_type, const Matrix<Integer>& Full_Cone_Generators);

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
