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

#include <cstdlib>
#include <csignal>
#include "libnormaliz/general.h"

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

int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (long) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);
    return 0;
}
#endif

namespace libnormaliz {

bool verbose = false;
bool constructor_verbose = true;
bool polynomial_verbose = false;
bool talkative = false;

volatile sig_atomic_t nmz_interrupted = 0;
const int default_thread_limit = 8;
int thread_limit = default_thread_limit;
bool parallelization_set = false;
bool no_output_on_interrupt = false;

// bool test_arithmetic_overflow = false;
// long overflow_test_modulus = 15401;

size_t GMP_mat = 0;
size_t GMP_hyp = 0;
size_t GMP_scal_prod = 0;
size_t TotDet = 0;

long cone_recursion_level = 0;
long full_cone_recursion_level = 0;

bool int_max_value_dual_long_computed = false;
bool int_max_value_dual_long_long_computed = false;
bool int_max_value_primary_long_computed = false;
bool int_max_value_primary_long_long_computed = false;

vector<vector<vector<long> > > CollectedAutoms(default_thread_limit);  // for use in nmz_nauty.cpp

long split_index_option = -1;
long split_index_rounds = -1;
long split_refinement = -1;
bool is_split_patching = false;

bool list_of_input_files = false;
long number_normaliz_instances= -1;

std::string global_project;

#ifdef NMZ_EXTENDED_TESTS
bool test_arith_overflow_full_cone = false;
bool test_arith_overflow_dual_mode = false;
bool test_arith_overflow_descent = false;
bool test_arith_overflow_proj_and_lift = false;
bool test_simplex_parallel = false;
bool test_small_pyramids = false;
bool test_large_pyramids = false;
bool test_linear_algebra_GMP = false;
#endif

#ifdef NMZ_NAUTY
void kill_nauty();
#endif

void interrupt_signal_handler(int signal) {
    nmz_interrupted = 1;
#ifdef NMZ_NAUTY
    kill_nauty();
#endif
}

namespace {
std::ostream* verbose_ostream_ptr = &std::cout;
std::ostream* error_ostream_ptr = &std::cerr;
}  // namespace

bool setVerboseDefault(bool v) {
    // we want to return the old value
    bool old = verbose;
    verbose = v;
    return old;
}

bool setTalkativeDefault(bool v) {
    // we want to return the old value
    bool old_talk = talkative;
    talkative = v;
    return old_talk;
}

bool setPolynomialVerbose(bool onoff){
    bool old = polynomial_verbose;
    polynomial_verbose = onoff;
    return old;
}

void suppressNextConstructorVerbose(){
        constructor_verbose = false;
}

int set_thread_limit(int t) {
    int old = thread_limit;
    parallelization_set = true;
    thread_limit = t;
    CollectedAutoms.resize(t);
    return old;
}

void setVerboseOutput(std::ostream& v_out) {
    verbose_ostream_ptr = &v_out;
}

void setErrorOutput(std::ostream& e_out) {
    error_ostream_ptr = &e_out;
}

std::ostream& verboseOutput() {
    return *verbose_ostream_ptr;
}

std::ostream& errorOutput() {
    return *error_ostream_ptr;
}

struct timeval TIME_global_begin, TIME_step_begin;
double GlobalTimeBound = -1.0; // can be set in normaliz.cfg
double GlobalPredictionTimeBound = -1.0;
// -1.0 means: no time bound

void StartTime(struct timeval& var_TIME_begin) {
    gettimeofday(&var_TIME_begin, 0);
}

void StartTime(){
        StartTime(TIME_step_begin);
}

void StartGlobalTime() {
    StartTime(TIME_global_begin);
}

double MeasureTime(const struct timeval var_TIME_begin) {
    struct timeval time_end;
    gettimeofday(&time_end, 0);
    long seconds = time_end.tv_sec - var_TIME_begin.tv_sec;
    long microseconds = time_end.tv_usec - var_TIME_begin.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    return elapsed;
}

double TimeSinceStart(){
    double elapsed = MeasureTime(TIME_global_begin);
    return elapsed;
}

void MeasureGlobalTime(bool verbose) {
    double elapsed = TimeSinceStart();
    if (verbose)
        verboseOutput() << "Normaliz elapsed wall clock time: " << elapsed << " sec" << endl;
}

void PrintTime(const struct timeval var_TIME_begin, bool verbose, const std::string& step){
    double elapsed = MeasureTime(var_TIME_begin);
    if (verbose)
        verboseOutput() << step << ": " << elapsed << " sec" << endl;
}

void MeasureTime(bool verbose, const std::string& step){
    PrintTime(TIME_step_begin,verbose,step);
}

unsigned int getVersion()
{
#ifdef NMZ_MAKEFILE_CLASSIC
#ifdef NMZ_DEVELOP
#define NMZ_RELEASE 9999999
#endif
#endif
    return NMZ_RELEASE;
}

} /* end namespace libnormaliz */
