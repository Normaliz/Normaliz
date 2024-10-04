/*
 * Normaliz
 * Copyright (C) 2007-2022  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#ifndef LIBNORMALIZ_GENERAL_H_
#define LIBNORMALIZ_GENERAL_H_

#include <iostream>
#include <cassert>
#include <csignal>
#include <cstddef>
#include <string>
#include <vector>

#ifndef _MSC_VER
#include <sys/time.h>
#else
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64

// MSVC defines this in winsock2.h!?
typedef struct timeval {
    long tv_sec;
    long tv_usec;
} timeval;
#endif

#include <libnormaliz/dynamic_bitset.h>

#ifndef NMZ_MAKEFILE_CLASSIC
#include <libnormaliz/nmz_config.h>
#endif

#ifdef _WIN32
#if defined(DLL_EXPORT)
#define NORMALIZ_DLL_EXPORT __declspec(dllexport)
#elif defined(NORMALIZ_USE_DLL) && !defined(NORMALIZ_USE_STATIC)
#define NORMALIZ_DLL_EXPORT __declspec(dllimport)
#else
#define NORMALIZ_DLL_EXPORT
#endif
#else
#define NORMALIZ_DLL_EXPORT
#endif

#ifndef NMZ_DEVELOP
#include "libnormaliz/version.h"
#endif

#include "libnormaliz/my_omp.h"

#ifdef USE_MPIR
#include <mpirxx.h>
#else  // otherwise use GMP
#include <gmpxx.h>
#endif

// in the serial version there is no need to catch-rethrow
#ifndef _OPENMP
#define NCATCH
#endif

#ifdef ENFNORMALIZ
#include <e-antic/renfxx.h>
namespace libnormaliz {
using eantic::renf_class;
using eantic::renf_elem_class;
typedef boost::intrusive_ptr<const renf_class> renf_class_shared;
}  // namespace libnormaliz
#else
namespace libnormaliz {
typedef long renf_elem_class;
struct renf_class {};
typedef renf_class* renf_class_shared;
}  // namespace libnormaliz
#endif

namespace libnormaliz {

typedef long long MachineInteger;
typedef double nmz_float;
const nmz_float nmz_epsilon = 1.0e-12;

/* this type is used in the entries of keys
 * it has to be able to hold number of generators */
typedef unsigned int key_t;
typedef unsigned short shortkey_t;

NORMALIZ_DLL_EXPORT extern bool verbose;
NORMALIZ_DLL_EXPORT extern bool running_input_file;
NORMALIZ_DLL_EXPORT extern bool constructor_verbose;
NORMALIZ_DLL_EXPORT extern bool polynomial_verbose;
NORMALIZ_DLL_EXPORT extern bool talkative;
NORMALIZ_DLL_EXPORT extern size_t GMP_mat, GMP_hyp, GMP_scal_prod;
NORMALIZ_DLL_EXPORT extern size_t TotDet;

NORMALIZ_DLL_EXPORT extern bool int_max_value_dual_long_computed;
NORMALIZ_DLL_EXPORT extern bool int_max_value_dual_long_long_computed;
NORMALIZ_DLL_EXPORT extern bool int_max_value_primary_long_computed;
NORMALIZ_DLL_EXPORT extern bool int_max_value_primary_long_long_computed;
NORMALIZ_DLL_EXPORT extern bool no_output_on_interrupt;
NORMALIZ_DLL_EXPORT extern bool no_coord_transf;
NORMALIZ_DLL_EXPORT extern bool write_lp_file;

#ifdef NMZ_EXTENDED_TESTS
NORMALIZ_DLL_EXPORT extern bool test_arith_overflow_full_cone, test_arith_overflow_dual_mode;
NORMALIZ_DLL_EXPORT extern bool test_arith_overflow_descent, test_arith_overflow_proj_and_lift;
NORMALIZ_DLL_EXPORT extern bool test_small_pyramids, test_large_pyramids;
NORMALIZ_DLL_EXPORT extern bool test_linear_algebra_GMP, test_simplex_parallel;
#endif

/*
 * If this variable is set to true, the current computation is interrupted and
 * an InterruptException is raised.
 */
NORMALIZ_DLL_EXPORT extern volatile sig_atomic_t nmz_interrupted;

// NORMALIZ_DLL_EXPORT extern bool nmz_scip; // controls the use of Scip

#define INTERRUPT_COMPUTATION_BY_EXCEPTION              \
    if (nmz_interrupted) {                              \
        throw InterruptException("external interrupt"); \
    }

/* if test_arithmetic_overflow is true, many operations are also done
 * modulo overflow_test_modulus to ensure the correctness of the calculations */
// extern bool test_arithmetic_overflow;
// extern long overflow_test_modulus;

NORMALIZ_DLL_EXPORT extern const int default_thread_limit;
NORMALIZ_DLL_EXPORT extern int thread_limit;
NORMALIZ_DLL_EXPORT extern bool parallelization_set;
int set_thread_limit(int t);

// debugging helpers
NORMALIZ_DLL_EXPORT extern long cone_recursion_level;
NORMALIZ_DLL_EXPORT extern long full_cone_recursion_level;

/* set the verbose default value */
bool setVerboseDefault(bool v);
void suppressNextConstructorVerbose();
bool setPolynomialVerbose(bool onoff);
bool setTalkativeDefault(bool v);
/* methods to set and use the output streams */
void setVerboseOutput(std::ostream&);
void setErrorOutput(std::ostream&);

std::ostream& verboseOutput();
std::ostream& errorOutput();

void interrupt_signal_handler(int signal);

void StartTime(struct timeval& var_TIME_begin);
void StartGlobalTime();
void StartTime();
double MeasureTime(const struct timeval var_TIME_begin);
void MeasureTime(bool verbose, const std::string& step);
double TimeSinceStart();
void MeasureGlobalTime(bool verbose);
void PrintTime(const struct timeval var_TIME_begin, bool verbose, const std::string& step);
void noCoordTransf(bool onoff);

void Check_Stop();


NORMALIZ_DLL_EXPORT extern double GlobalTimeBound;
NORMALIZ_DLL_EXPORT extern double GlobalPredictionTimeBound;


NORMALIZ_DLL_EXPORT extern long level_local_solutions; // transports <l> of -Q=<n>
NORMALIZ_DLL_EXPORT extern long split_index_option; // transports <n> of -X=<n>
NORMALIZ_DLL_EXPORT extern long split_index_rounds; // transports the split index option after adding the rounds
NORMALIZ_DLL_EXPORT extern long split_refinement; // transports the refinement of the split
NORMALIZ_DLL_EXPORT extern bool is_split_patching; // indicates that we are computing a split
NORMALIZ_DLL_EXPORT extern bool save_local_solutions; // indicates that local solutions are to be stored in distributed computation

NORMALIZ_DLL_EXPORT extern bool list_of_input_files; // true if processing list of input files
NORMALIZ_DLL_EXPORT extern long number_normaliz_instances; // for distribution of input files to several instances of normaliz
NORMALIZ_DLL_EXPORT extern long input_file_option; // index modulo number_normaliz_instances  of input fileas to be run by this instance
NORMALIZ_DLL_EXPORT extern size_t verb_length; // helps in verrbose output

NORMALIZ_DLL_EXPORT extern std::string global_project;
NORMALIZ_DLL_EXPORT extern std::string lat_file_name;

// The following hold data read from the input file, used inside fusion.cp
/*
NORMALIZ_DLL_EXPORT extern std::vector<key_t> fusion_type_coinc_from_input;
NORMALIZ_DLL_EXPORT extern std::string fusion_type_from_input;
NORMALIZ_DLL_EXPORT extern std::vector<key_t> fusion_duality_from_input;
NORMALIZ_DLL_EXPORT extern std::vector<key_t> candidate_subring_from_input;
NORMALIZ_DLL_EXPORT extern std::vector<key_t> fusion_type_for_partition_from_input;
NORMALIZ_DLL_EXPORT extern bool fusion_commutative_from_input;
*/
NORMALIZ_DLL_EXPORT extern bool write_fusion_mult_tables_from_input;
// void set_global_fusion_data();
// void reset_global_fusion_data();
/*
NORMALIZ_DLL_EXPORT extern bool write_fusion_mult_tables_from_input;
void set_global_fusion_data();
void reset_global_fusion_data();
*/


} /* end namespace libnormaliz */

#include <libnormaliz/normaliz_exception.h>
// #include <libnormaliz/input_type.h>
#include <libnormaliz/cone_property.h>
#include <libnormaliz/integer.h>

#endif /* GENERAL_H_ */
