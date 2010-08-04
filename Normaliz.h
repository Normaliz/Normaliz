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
#ifndef NORMALIZ_H
#define NORMALIZ_H
//---------------------------------------------------------------------------

#include <assert.h>

#ifndef NO_OPENMP
#include <omp.h>
#endif

#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "simplex.h"
#include "list_operations.h"
#include "lineare_transformation.h"
#include "full_cone.h"
#include "output.h"
#include "mode.h"



//---------------------------------------------------------------------------
// global variables

// used for turn on and off the tests for arithmetic overflow
// the run time may double when the tests are performed
bool test_arithmetic_overflow=false;

// used for testing possible arithmetic overflow at key points
// the bigger the number test is, the bigger the probability that no arithmetic
// overflow is produced
// a matrix which procuses overflow over int but not over long long is
//    22222   33337
//    55559   77773
int overflow_test_modulus=15401;

// used to turn on and off the display of a progress report
// designed for users who run complex examples
// it is activated by setup or by the option 'c' in the command line
bool verbose=false;

// used to turn on and off speed optimization
// depend on how much RAM is available on the system
bool optimize_speed=true;


//---------------------------------------------------------------------------
// global functions

/**
 * Determinates if and how the program will be terminated in case of errors
 */
void global_error_handling();

/**
 * Prints help text
 * @param command Name of the executable
 */
void printHelp(char* command);

int main(int argc, char* argv[]);

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
