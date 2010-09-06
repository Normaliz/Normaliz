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

#include "libnormaliz.h"

#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "list_operations.h"
#include "lineare_transformation.h"
#include "full_cone.h"
#include "output.h"
#include "mode.h"


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
