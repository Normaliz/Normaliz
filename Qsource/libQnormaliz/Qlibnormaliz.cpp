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
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include "libQnormaliz/libQnormaliz.h"
#include "libQnormaliz/Qgeneral.h"

namespace libQnormaliz {

bool verbose = false;

volatile sig_atomic_t nmz_interrupted = 0;
long default_thread_limit=8;
long thread_limit=default_thread_limit;
bool parallelization_set=false;

long set_thread_limit(long t){
    long old=thread_limit;
    parallelization_set=true;
    thread_limit=t;
    return old;
}

// bool test_arithmetic_overflow = false;
// long overflow_test_modulus = 15401;

size_t GMP_mat=0;
size_t GMP_hyp=0;
size_t GMP_scal_prod=0;
size_t TotDet=0;

void interrupt_signal_handler( int signal ){
    nmz_interrupted = 1;
}


namespace {
    std::ostream* verbose_ostream_ptr = &std::cout;
    std::ostream* error_ostream_ptr = &std::cerr;
} // end anonymous namespace, only accessible in this file (and when it is included)

bool setVerboseDefault(bool v) {
    //we want to return the old value
    bool old = verbose;
    verbose = v;
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

InputType to_type(const std::string& type_string) {

    if ( type_string=="0" || type_string=="1" || type_string=="2" || type_string=="3"
      || type_string=="4" || type_string=="5" || type_string=="6"
      || type_string=="hyperplanes"
      || type_string=="10") {
        throw BadInputException("Error: deprecated type \"" + type_string
                + "\", please use new type string!");
    }

    if (type_string=="0"||type_string=="integral_closure") {
        return QType::integral_closure;
    }
    if (type_string=="polyhedron") {
        return QType::polyhedron;
    }
    if (type_string=="1"||type_string=="normalization") {
        return QType::normalization;
    }
    if (type_string=="2"||type_string=="polytope") {
        return QType::polytope;
    }
    if (type_string=="3"||type_string=="rees_algebra") {
        return QType::rees_algebra;
    }
    if (type_string=="4"||type_string=="hyperplanes" ||type_string=="inequalities") {
        return QType::inequalities;
    }
    if (type_string=="strict_inequalities") {
         return QType::strict_inequalities;
    }
    if (type_string=="strict_signs") {
        return QType::strict_signs;
    }
    if (type_string=="inhom_inequalities") {
        return QType::inhom_inequalities;
    }
    if (type_string=="dehomogenization") {
         return QType::dehomogenization;
    }
    if (type_string=="5"||type_string=="equations") {
        return QType::equations;
    }
    if (type_string=="inhom_equations") {
        return QType::inhom_equations;
    }
    if (type_string=="6"||type_string=="congruences") {
        return QType::congruences;
    }
    if (type_string=="inhom_congruences") {
        return QType::inhom_congruences;
    }
    if (type_string=="signs") {
        return QType::signs;
    }
    if (type_string=="10"||type_string=="lattice_ideal") {
        return QType::lattice_ideal;
    }
    if (type_string=="grading") {
        return QType::grading;
    }
    if (type_string=="excluded_faces") {
        return QType::excluded_faces;
    }
    if (type_string=="lattice") {
        return QType::lattice;
    }
    if (type_string=="saturation") {
        return QType::saturation;
    }
    if (type_string=="cone") {
        return QType::cone;
    }
    if (type_string=="offset") {
        return QType::offset;
    }
    if (type_string=="vertices") {
        return QType::vertices;
    }
    if (type_string=="support_hyperplanes") {
        return QType::support_hyperplanes;
    }
    if (type_string=="cone_and_lattice") {
        return QType::cone_and_lattice;
    }
    if (type_string=="subspace") {
        return QType::subspace;
    }

    throw BadInputException("Unknown type \"" + type_string + "\"!");
    return QType::integral_closure;
}

long type_nr_columns_correction(InputType t) {
    if (t == QType::polytope || t == QType::rees_algebra)
        return -1;
    if (t == QType::congruences || t == QType::vertices || t == QType::polyhedron
     || t == QType::inhom_inequalities || t == QType::inhom_equations)
        return 1;
    if (t == QType::inhom_congruences)
        return 2;
    return 0;
}

/* returns true if the input of this type is a vector */
bool type_is_vector(InputType type){
    if (type == QType::grading || type == QType::signs || type == QType::strict_signs
            || type == QType::dehomogenization || type == QType::offset) {
        return true;
    }
    return false;
}

} /* end namespace libQnormaliz */
