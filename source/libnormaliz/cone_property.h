/*
 * Normaliz
 * Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#ifndef LIBNORMALIZ_CONE_PROPERTY_H_
#define LIBNORMALIZ_CONE_PROPERTY_H_

#include <bitset>
#include <ostream>

namespace libnormaliz {

/****************************************************************************
**
**  'START_ENUM_RANGE' and 'END_ENUM_RANGE' simplify creating "ranges" of
**  enum variables.
**
**  Usage example:
**    enum {
**      START_ENUM_RANGE(FIRST),
**        FOO,
**        BAR,
**      END_ENUM_RANGE(LAST)
**    };
**  is essentially equivalent to
**    enum {
**      FIRST,
**        FOO = FIRST,
**        BAR,
**      LAST = BAR
**    };
**  Note that if we add a value into the range after 'BAR', we must adjust
**  the definition of 'LAST', which is easy to forget. Also, reordering enum
**  values may require extra work. With the range macros, all of this is
**  taken care of automatically.
*/
#define START_ENUM_RANGE(id)            id, _##id##_post = id - 1
#define END_ENUM_RANGE(id)              _##id##_pre, id = _##id##_pre - 1


/* An enumeration of things, that can be computed for a cone.
 * The namespace prevents interfering with other names.
 * Remember to change also the string conversion if you change this enum!
 */
namespace ConeProperty {
    enum Enum {
        // matrix valued
        START_ENUM_RANGE(FIRST_MATRIX),
        Generators,
        ExtremeRays,
        VerticesOfPolyhedron,
        SupportHyperplanes,
        HilbertBasis,
        ModuleGenerators,
        Deg1Elements,
        LatticePoints,
        ModuleGeneratorsOverOriginalMonoid,
        ExcludedFaces,
        OriginalMonoidGenerators,
        MaximalSubspace,
        Equations,
        Congruences,
        END_ENUM_RANGE(LAST_MATRIX),

        START_ENUM_RANGE(FIRST_MATRIX_FLOAT),
        SuppHypsFloat,
        VerticesFloat,
        END_ENUM_RANGE(LAST_MATRIX_FLOAT),

        // vector valued
        START_ENUM_RANGE(FIRST_VECTOR),
        Grading,
        Dehomogenization,
        WitnessNotIntegrallyClosed,
        GeneratorOfInterior,
        ClassGroup,
        END_ENUM_RANGE(LAST_VECTOR),

        // integer valued
        START_ENUM_RANGE(FIRST_INTEGER),
        TriangulationDetSum,
        ReesPrimaryMultiplicity,
        GradingDenom,
        UnitGroupIndex,
        InternalIndex,
        END_ENUM_RANGE(LAST_INTEGER),

        START_ENUM_RANGE(FIRST_GMP_INTEGER),
        ExternalIndex = FIRST_GMP_INTEGER,
        END_ENUM_RANGE(LAST_GMP_INTEGER),

        // rational valued
        START_ENUM_RANGE(FIRST_RATIONAL),
        Multiplicity,
        Volume,
        Integral,
        VirtualMultiplicity,
        END_ENUM_RANGE(LAST_RATIONAL),

        // field valued
        START_ENUM_RANGE(FIRST_FIELD_ELEM),
        RenfVolume=FIRST_FIELD_ELEM,
        LAST_FIELD_ELEM=ConeProperty::RenfVolume,

        // floating point valued
        START_ENUM_RANGE(FIRST_FLOAT),
        EuclideanVolume,
        EuclideanIntegral,
        END_ENUM_RANGE(LAST_FLOAT),

        // dimensions
        START_ENUM_RANGE(FIRST_MACHINE_INTEGER),
        TriangulationSize,
        NumberLatticePoints,
        RecessionRank,
        AffineDim,
        ModuleRank,
        Rank,
        EmbeddingDim,
        END_ENUM_RANGE(LAST_MACHINE_INTEGER),

        // boolean valued 
        START_ENUM_RANGE(FIRST_BOOLEAN),
        IsPointed,
        IsDeg1ExtremeRays,
        IsDeg1HilbertBasis,
        IsIntegrallyClosed,
        IsReesPrimary,
        IsInhomogeneous,
        IsGorenstein,
        END_ENUM_RANGE(LAST_BOOLEAN),

        // complex structures
        START_ENUM_RANGE(FIRST_COMPLEX_STRUCTURE),
        Triangulation,
        StanleyDec,
        InclusionExclusionData,
        IntegerHull,
        ProjectCone,
        ConeDecomposition,
        //
        Automorphisms,
        AmbientAutomorphisms,
        CombinatorialAutomorphisms,
        RationalAutomorphisms,
        EuclideanAutomorphisms,
        //
        HilbertSeries,
        HilbertQuasiPolynomial,
        EhrhartSeries,
        EhrhartQuasiPolynomial,
        WeightedEhrhartSeries,
        WeightedEhrhartQuasiPolynomial,
        FaceLattice,
        FVector,
        Incidence,
        Sublattice,
        END_ENUM_RANGE(LAST_COMPLEX_STRUCTURE),

        //
        // integer type for computations
        //
        START_ENUM_RANGE(FIRST_PROPERTY),
        BigInt,
        //
        // algorithmic variants
        //
        DefaultMode,
        Approximate,
        BottomDecomposition,
        NoBottomDec,       
        DualMode,
        PrimalMode,
        Projection,
        ProjectionFloat,
        NoProjection,
        Symmetrize,
        NoSymmetrization,
        NoSubdivision,
        NoNestedTri, // synonym for NoSubdivision
        KeepOrder,
        HSOP,
        NoPeriodBound,
        NoLLL,
        NoRelax,
        Descent,
        NoDescent,
        NoGradingDenom,
        GradingIsPositive,
        ExploitAutomsVectors,
        ExploitAutomsMult,
        //
        Dynamic,
        Static,
        //
        // checking properties of already computed data
        // (cannot be used as a computation goal)
        //
        IsTriangulationNested,
        IsTriangulationPartial,
        //
        // ONLY FOR INTERNAL CONTROL
        //
        // ExplicitHilbertSeries,
        NakedDual,
        END_ENUM_RANGE(LAST_PROPERTY),

        EnumSize // this has to be the last entry, to get the number of entries in the enum

    }; // remember to change also the string conversion function if you change this enum
}

namespace OutputType{
    enum Enum {
        Matrix,
        MatrixFloat,
        Vector,
        Integer,
        GMPInteger,
        Rational,
        FieldElem,
        Float,
        MachineInteger,
        Bool,
        Complex,
        Void
    };
}

class ConeProperties {
public:
    /* Constructors */
    ConeProperties();
    ConeProperties(ConeProperty::Enum);
    ConeProperties(ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties(ConeProperty::Enum, ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties(const std::bitset<ConeProperty::EnumSize>&);

    /* set properties */
    ConeProperties& set(ConeProperty::Enum, bool value=true);
    ConeProperties& set(const std::string s, bool value=true);
    ConeProperties& set(ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties& set(ConeProperty::Enum, ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties& set(const ConeProperties&);

    /* reset (=unset) properties */
    ConeProperties& reset(ConeProperty::Enum Property);
    ConeProperties& reset(const ConeProperties&);
    ConeProperties& reset_compute_options();

    /* test which/how many properties are set */
    bool test(ConeProperty::Enum Property) const;
    bool any() const;
    bool none() const;
    size_t count () const;

    /* return the restriction of this to the goals / options */
    ConeProperties goals();
    ConeProperties options();

    /* the following methods are used internally */
    void set_preconditions(bool inhomogeneous, bool numberfield);    // activate properties which are needed implicitly
    // void prepare_compute_options(bool inhomogeneous, bool numberfield);
    void check_sanity(bool inhomogeneous);

    void check_conflicting_variants();
    void check_Q_permissible(bool after_implications);
    // void set_default_goals(bool inhomogeneous, bool numberfield);

    /* print it in a nice way */
    friend std::ostream& operator<<(std::ostream&, const ConeProperties&);


private:
    std::bitset<ConeProperty::EnumSize> CPs;

};

// conversion to/from strings
bool isConeProperty(ConeProperty::Enum& cp, const std::string& s);
ConeProperty::Enum toConeProperty(const std::string&);
const std::string& toString(ConeProperty::Enum);
std::ostream& operator<<(std::ostream&, const ConeProperties&);
OutputType::Enum output_type(ConeProperty::Enum);

}

#endif /* CONE_PROPERTY_H_ */
