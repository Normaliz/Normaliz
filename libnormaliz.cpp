/*
 * libnormaliz.cpp
 *
 *  Created on: Aug 24, 2010
 *      Author: csoeger
 */

#include "libnormaliz.h"



bool verbose = false;

bool test_arithmetic_overflow = false;
int overflow_test_modulus = 15401;

std::ostream& verbose_ostream = std::cout;
std::ostream& error_ostream = std::cerr;
