/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "libQnormaliz/libQnormaliz.h"
#include "libQnormaliz/Qcone_property.h"


void printHeader();
void printCopying();
void printVersion();

/**
 * Prints help text
 * @param command Name of the executable
 */
void printHelp(char* command);

int main(int argc, char* argv[]);

/*template<typename Integer>
int process_data(string& output_name, libnormaliz::ConeProperties to_compute, bool write_extra_files, bool write_all_files, bool verbose);
*/

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
