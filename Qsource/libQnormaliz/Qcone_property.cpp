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

#ifdef NMZ_MIC_OFFLOAD
#pragma offload_attribute (push, target(mic))
#endif

#include <vector>
#include <string>
#include <assert.h>

#include "libQnormaliz/Qcone_property.h"
#include "libQnormaliz/libQnormaliz.h"
#include "libQnormaliz/Qnormaliz_exception.h"

namespace libQnormaliz {
using std::bitset;
using std::vector;
using std::string;
using std::endl;


/* Constructors */
ConeProperties::ConeProperties() {
    CPs = bitset<QConeProperty::EnumSize>();
}
ConeProperties::ConeProperties(QConeProperty::Enum p1) {
    CPs = bitset<QConeProperty::EnumSize>();
    CPs.set(p1);
}
ConeProperties::ConeProperties(QConeProperty::Enum p1, QConeProperty::Enum p2) {
    CPs = bitset<QConeProperty::EnumSize>();
    CPs.set(p1);
    CPs.set(p2);
}
ConeProperties::ConeProperties(QConeProperty::Enum p1, QConeProperty::Enum p2,
                               QConeProperty::Enum p3) {
    CPs = bitset<QConeProperty::EnumSize>();
    CPs.set(p1);
    CPs.set(p2);
    CPs.set(p3);
}
ConeProperties::ConeProperties(const bitset<QConeProperty::EnumSize>& props){
    CPs = props;
}

/* set Properties */
ConeProperties& ConeProperties::set(QConeProperty::Enum p1, bool value) {
    CPs.set(p1, value);
    return *this;
}
ConeProperties& ConeProperties::set(QConeProperty::Enum p1, QConeProperty::Enum p2) {
    CPs.set(p1);
    CPs.set(p2);
    return *this;
}
ConeProperties& ConeProperties::set(const ConeProperties& ConeProps) {
    CPs ^= ConeProps.CPs;
    return *this;
}

ConeProperties& ConeProperties::set(const std::string s, bool value) {
    CPs.set(toConeProperty(s), value);
    return *this;
}

/* reset (=unset) properties */
ConeProperties& ConeProperties::reset(QConeProperty::Enum Property) {
    CPs.set(Property, false);
    return *this;
}
ConeProperties& ConeProperties::reset(const ConeProperties& ConeProps) {
    CPs &= ~ConeProps.CPs;
    return *this;
}

ConeProperties& ConeProperties::reset_compute_options() {
    CPs.set(QConeProperty::Projection, false);
    CPs.set(QConeProperty::ProjectionFloat, false);
    CPs.set(QConeProperty::NoProjection, false);
    CPs.set(QConeProperty::Approximate, false);
    CPs.set(QConeProperty::BottomDecomposition, false);
    CPs.set(QConeProperty::NoBottomDec, false);
    CPs.set(QConeProperty::DefaultMode, false);
    CPs.set(QConeProperty::DualMode, false);
    CPs.set(QConeProperty::PrimalMode, false);
    CPs.set(QConeProperty::KeepOrder, false);
    CPs.set(QConeProperty::HSOP, false);
    CPs.set(QConeProperty::Symmetrize, false);
    CPs.set(QConeProperty::NoSymmetrization, false);
    CPs.set(QConeProperty::BigInt, false);
    CPs.set(QConeProperty::NoSubdivision, false);
    CPs.set(QConeProperty::NoNestedTri, false);
    CPs.set(QConeProperty::NoPeriodBound, false);
    CPs.set(QConeProperty::SCIP, false);
    CPs.set(QConeProperty::NoLLL, false);
    CPs.set(QConeProperty::NoRelax, false);
    CPs.set(QConeProperty::ExplicitHilbertSeries, false);
    CPs.set(QConeProperty::NakedDual, false);
    CPs.set(QConeProperty::Descent, false);
    CPs.set(QConeProperty::NoDescent, false);
    CPs.set(QConeProperty::NoGradingDenom, false);
    CPs.set(QConeProperty::GradingIsPositive, false);
    return *this;
}


/* return a new ConeProperties object with only the goals/options set,
 * which are set in this object
 */
ConeProperties ConeProperties::goals() {
    ConeProperties ret(*this);
    ret.reset_compute_options();
    return ret;
}

ConeProperties ConeProperties::options() {
    ConeProperties ret;
    ret.set(QConeProperty::Projection, CPs.test(QConeProperty::Projection));
    ret.set(QConeProperty::ProjectionFloat, CPs.test(QConeProperty::ProjectionFloat));
    ret.set(QConeProperty::NoProjection, CPs.test(QConeProperty::NoProjection));
    ret.set(QConeProperty::Approximate, CPs.test(QConeProperty::Approximate));
    ret.set(QConeProperty::BottomDecomposition, CPs.test(QConeProperty::BottomDecomposition));
    ret.set(QConeProperty::NoBottomDec, CPs.test(QConeProperty::NoBottomDec));
    ret.set(QConeProperty::DefaultMode, CPs.test(QConeProperty::DefaultMode));
    ret.set(QConeProperty::DualMode, CPs.test(QConeProperty::DualMode));
    ret.set(QConeProperty::KeepOrder, CPs.test(QConeProperty::KeepOrder));
    ret.set(QConeProperty::HSOP, CPs.test(QConeProperty::HSOP));
    ret.set(QConeProperty::Symmetrize, CPs.test(QConeProperty::Symmetrize));
    ret.set(QConeProperty::NoSymmetrization, CPs.test(QConeProperty::NoSymmetrization));
    ret.set(QConeProperty::PrimalMode, CPs.test(QConeProperty::PrimalMode));
    ret.set(QConeProperty::NoSubdivision, CPs.test(QConeProperty::NoSubdivision));
    ret.set(QConeProperty::NoNestedTri, CPs.test(QConeProperty::NoNestedTri));
    ret.set(QConeProperty::BigInt, CPs.test(QConeProperty::BigInt));
    ret.set(QConeProperty::NoPeriodBound, CPs.test(QConeProperty::NoPeriodBound));
    ret.set(QConeProperty::SCIP, CPs.test(QConeProperty::SCIP));
    ret.set(QConeProperty::NoLLL, CPs.test(QConeProperty::NoLLL));
    ret.set(QConeProperty::NoRelax, CPs.test(QConeProperty::NoRelax));
    ret.set(QConeProperty::ExplicitHilbertSeries, CPs.test(QConeProperty::ExplicitHilbertSeries));
    ret.set(QConeProperty::NakedDual, CPs.test(QConeProperty::NakedDual));
    ret.set(QConeProperty::Descent, CPs.test(QConeProperty::Descent));
    ret.set(QConeProperty::NoDescent, CPs.test(QConeProperty::NoDescent));
    ret.set(QConeProperty::NoGradingDenom, CPs.test(QConeProperty::NoGradingDenom));
    ret.set(QConeProperty::GradingIsPositive, CPs.test(QConeProperty::GradingIsPositive));
    return ret;
}

/* test which/how many properties are set */
bool ConeProperties::test(QConeProperty::Enum Property) const {
    return CPs.test(Property);
}
bool ConeProperties::any() const {
    return CPs.any();
}
bool ConeProperties::none() const {
    return CPs.none();
}
size_t ConeProperties::count() const {
    return CPs.count();
}


/* add preconditions */
void ConeProperties::set_preconditions(bool inhomogeneous) {
    
    if(inhomogeneous && CPs.test(QConeProperty::Deg1Elements)){
        CPs.set(QConeProperty::ModuleGenerators);
        CPs.reset(QConeProperty::Deg1Elements);
    }
    
    if(inhomogeneous && CPs.test(QConeProperty::LatticePoints)){
        //CPs.set(QConeProperty::ModuleGenerators);
        CPs.set(QConeProperty::HilbertBasis);
        CPs.reset(QConeProperty::Deg1Elements);
        CPs.reset(QConeProperty::LatticePoints);
    }
    
    if(!inhomogeneous &&  CPs.test(QConeProperty::LatticePoints)){
        CPs.set(QConeProperty::NoGradingDenom);
        CPs.set(QConeProperty::Deg1Elements);
        CPs.reset(QConeProperty::LatticePoints);
    }
    
    if(CPs.test(QConeProperty::EuclideanVolume))
        CPs.set(QConeProperty::Volume);

    if (CPs.test(QConeProperty::IsDeg1ExtremeRays)) {
        CPs.set(QConeProperty::ExtremeRays);
        CPs.set(QConeProperty::Grading);
    }
    if (CPs.test(QConeProperty::Grading))
        CPs.set(QConeProperty::Generators);

    if (CPs.test(QConeProperty::IsPointed))
        CPs.set(QConeProperty::ExtremeRays);

    if (CPs.test(QConeProperty::ExtremeRays))
        CPs.set(QConeProperty::SupportHyperplanes);

    // inhomogenous preconditions
    if (CPs.test(QConeProperty::VerticesOfPolyhedron))
        CPs.set(QConeProperty::ExtremeRays);
    
    if (CPs.test(QConeProperty::MaximalSubspace))
        CPs.set(QConeProperty::SupportHyperplanes);
    
    // always
    
    if (CPs.test(QConeProperty::ExtremeRays))
        CPs.set(QConeProperty::SupportHyperplanes);
}

/* removes ignored compute options and sets implications */
void ConeProperties::prepare_compute_options(bool inhomogeneous) {
    if (CPs.test(QConeProperty::IntegerHull)){
        if(inhomogeneous){
            CPs.set(QConeProperty::ModuleGenerators);
        }
        else{
            CPs.set(QConeProperty::Deg1Elements);
        }
    }
    // -d without -1 means: compute Hilbert basis in dual mode
    if (CPs.test(QConeProperty::DualMode) && !CPs.test(QConeProperty::Deg1Elements)){
        CPs.set(QConeProperty::HilbertBasis);
    }
    
    if(CPs.test(QConeProperty::ModuleGeneratorsOverOriginalMonoid)) // can't be computed in dual mode
        CPs.reset(QConeProperty::DualMode);

    // dual mode has priority, approximation makes no sense if HB is computed
    if(CPs.test(QConeProperty::DualMode) || CPs.test(QConeProperty::HilbertBasis))
        CPs.reset(QConeProperty::Approximate);

    if ((CPs.test(QConeProperty::DualMode) || CPs.test(QConeProperty::Approximate))
        && (CPs.test(QConeProperty::HilbertSeries) || CPs.test(QConeProperty::StanleyDec))
         && !CPs.test(QConeProperty::HilbertBasis)){
        CPs.reset(QConeProperty::DualMode); //it makes no sense to compute only deg 1 elements in dual mode
        CPs.reset(QConeProperty::Approximate); // or by approximation if the
    }                                            // Stanley decomposition must be computed anyway
    if (CPs.test(QConeProperty::Approximate)
            && !CPs.test(QConeProperty::Deg1Elements)) {
        errorOutput() << "WARNING: Approximate is ignored since Deg1Elements is not set."<< std::endl;
    }
    if (CPs.test(QConeProperty::ConeDecomposition))
        CPs.set(QConeProperty::Triangulation); 
    
    if (CPs.test(QConeProperty::GradingDenom))
        CPs.reset(QConeProperty::Grading);
    
    if(CPs.test(QConeProperty::UnitGroupIndex))
        CPs.set(QConeProperty::HilbertBasis);
    
    if(CPs.test(QConeProperty::Equations) || CPs.test(QConeProperty::Congruences) || CPs.test(QConeProperty::ExternalIndex))
        CPs.set(QConeProperty::Sublattice);
    
    if(CPs.test(QConeProperty::Rank))
        CPs.set(QConeProperty::Sublattice);
    
    if(CPs.test(QConeProperty::HilbertQuasiPolynomial))
        CPs.set(QConeProperty::HilbertSeries);
    
    if(inhomogeneous && CPs.test(QConeProperty::SupportHyperplanes))
        CPs.set(QConeProperty::AffineDim);

    if(CPs.test(QConeProperty::DefaultMode)){
        /* CPs.set(QConeProperty::HilbertBasis);
        CPs.set(QConeProperty::HilbertSeries);
        if(!inhomogeneous)
            CPs.set(QConeProperty::ClassGroup);*/
        CPs.set(QConeProperty::SupportHyperplanes);        
    }
}

void ConeProperties::check_Q_permissible() {
    ConeProperties copy(*this);
    copy.reset(QConeProperty::SupportHyperplanes);
    copy.reset(QConeProperty::ExtremeRays);
    copy.reset(QConeProperty::VerticesOfPolyhedron);
    copy.reset(QConeProperty::KeepOrder);
    copy.reset(QConeProperty::Triangulation); 
    copy.reset(QConeProperty::ConeDecomposition);
    copy.reset(QConeProperty::DefaultMode);
    copy.reset(QConeProperty::Generators);
    copy.reset(QConeProperty::Sublattice);
    copy.reset(QConeProperty::MaximalSubspace);
    copy.reset(QConeProperty::Equations);
    copy.reset(QConeProperty::Dehomogenization);
    copy.reset(QConeProperty::Rank);
    copy.reset(QConeProperty::EmbeddingDim);
    copy.reset(QConeProperty::IsPointed);
    copy.reset(QConeProperty::IsInhomogeneous);
    copy.reset(QConeProperty::AffineDim);
    copy.reset(QConeProperty::ModuleGenerators);
    copy.reset(QConeProperty::Deg1Elements);
    copy.reset(QConeProperty::Volume);
    copy.reset(QConeProperty::IntegerHull);
    copy.reset(QConeProperty::Generators);
    copy.reset(QConeProperty::TriangulationDetSum);
    copy.reset(QConeProperty::LatticePoints);
//     copy.reset(QConeProperty::TriangulationSize);
    
    //bvverboseOutput() << copy << endl;
    if(copy.any()){
        verboseOutput() << copy << endl;
        throw BadInputException("Cone Property not allowed for field coefficients");
    }
}


void ConeProperties::check_sanity(bool inhomogeneous) {
    QConeProperty::Enum prop;
    if(        
           (CPs.test(QConeProperty::BottomDecomposition) && CPs.test(QConeProperty::NoBottomDec))
        || (CPs.test(QConeProperty::DualMode) && CPs.test(QConeProperty::PrimalMode))
        || (CPs.test(QConeProperty::Symmetrize) && CPs.test(QConeProperty::NoSymmetrization))
    )
        throw BadInputException("Contradictory algorithmic variants in options.");
        
    if(CPs.test(QConeProperty::IsTriangulationNested) || CPs.test(QConeProperty::IsTriangulationPartial))
        throw BadInputException("ConeProperty not allowed in compute().");
        
    for (size_t i=0; i<QConeProperty::EnumSize; i++) {
        if (CPs.test(i)) {
            prop = static_cast<QConeProperty::Enum>(i);
            if (inhomogeneous) {
                if ( prop == QConeProperty::Deg1Elements
                  || prop == QConeProperty::StanleyDec
                  // || prop == QConeProperty::Triangulation
                  || prop == QConeProperty::ConeDecomposition
                  || prop == QConeProperty::IsIntegrallyClosed
                  || prop == QConeProperty::WitnessNotIntegrallyClosed
                  || prop == QConeProperty::Approximate
                  || prop == QConeProperty::ClassGroup
                  || prop == QConeProperty::Symmetrize
                  || prop == QConeProperty::NoSymmetrization
                  || prop == QConeProperty::InclusionExclusionData
                  || prop == QConeProperty::ExcludedFaces
                  || prop == QConeProperty::UnitGroupIndex
                  || prop == QConeProperty::ReesPrimaryMultiplicity
                  || prop == QConeProperty::IsReesPrimary
                  || prop == QConeProperty::IsDeg1HilbertBasis
                  || prop == QConeProperty::IsDeg1ExtremeRays
                 // || prop == QConeProperty::ModuleGeneratorsOverOriginalMonoid
                ) {
                    throw BadInputException(toString(prop) + " not computable in the inhomogeneous case.");
                }
            } else { // homgeneous
                if ( prop == QConeProperty::VerticesOfPolyhedron
                  || prop == QConeProperty::ModuleRank
                  || prop == QConeProperty::ModuleGenerators ) {
                    throw BadInputException(toString(prop) + " only computable in the inhomogeneous case.");
                }
            }
        }  //end if test(i)
    }
}

/* conversion */
namespace {
    // only to initialize the CPN in ConePropertyNames
    vector<string> initializeCPN() {
        vector<string> CPN(QConeProperty::EnumSize);
        CPN.at(QConeProperty::Generators) = "Generators";
        CPN.at(QConeProperty::ExtremeRays) = "ExtremeRays";
        CPN.at(QConeProperty::VerticesFloat) = "VerticesFloat";
        CPN.at(QConeProperty::VerticesOfPolyhedron) = "VerticesOfPolyhedron";
        CPN.at(QConeProperty::SupportHyperplanes) = "SupportHyperplanes";
        CPN.at(QConeProperty::SuppHypsFloat) = "SuppHypsFloat";
        CPN.at(QConeProperty::TriangulationSize) = "TriangulationSize";
        CPN.at(QConeProperty::TriangulationDetSum) = "TriangulationDetSum";
        CPN.at(QConeProperty::Triangulation) = "Triangulation";
        CPN.at(QConeProperty::Multiplicity) = "Multiplicity";
        CPN.at(QConeProperty::Volume) = "Volume";
        CPN.at(QConeProperty::EuclideanVolume) = "EuclideanVolume";
        CPN.at(QConeProperty::RecessionRank) = "RecessionRank";
        CPN.at(QConeProperty::AffineDim) = "AffineDim";
        CPN.at(QConeProperty::ModuleRank) = "ModuleRank";
        CPN.at(QConeProperty::HilbertBasis) = "HilbertBasis";
        CPN.at(QConeProperty::ModuleGenerators) = "ModuleGenerators";
        CPN.at(QConeProperty::LatticePoints) = "LatticePoints";
        CPN.at(QConeProperty::Deg1Elements) = "Deg1Elements";
        CPN.at(QConeProperty::HilbertSeries) = "HilbertSeries";
        CPN.at(QConeProperty::Grading) = "Grading";
        CPN.at(QConeProperty::IsPointed) = "IsPointed";
        CPN.at(QConeProperty::IsDeg1ExtremeRays) = "IsDeg1ExtremeRays";
        CPN.at(QConeProperty::IsDeg1HilbertBasis) = "IsDeg1HilbertBasis";
        CPN.at(QConeProperty::IsIntegrallyClosed) = "IsIntegrallyClosed";
        CPN.at(QConeProperty::WitnessNotIntegrallyClosed) = "WitnessNotIntegrallyClosed";
        CPN.at(QConeProperty::OriginalMonoidGenerators) = "OriginalMonoidGenerators";
        CPN.at(QConeProperty::IsReesPrimary) = "IsReesPrimary";
        CPN.at(QConeProperty::ReesPrimaryMultiplicity) = "ReesPrimaryMultiplicity";
        CPN.at(QConeProperty::StanleyDec) = "StanleyDec";
        CPN.at(QConeProperty::ExcludedFaces) = "ExcludedFaces";
        CPN.at(QConeProperty::Dehomogenization) = "Dehomogenization";
        CPN.at(QConeProperty::InclusionExclusionData) = "InclusionExclusionData";
        CPN.at(QConeProperty::Sublattice) = "Sublattice";
        CPN.at(QConeProperty::ClassGroup) = "ClassGroup";
        CPN.at(QConeProperty::ModuleGeneratorsOverOriginalMonoid) = "ModuleGeneratorsOverOriginalMonoid";
        CPN.at(QConeProperty::Approximate) = "Approximate";
        CPN.at(QConeProperty::BottomDecomposition) = "BottomDecomposition";
        CPN.at(QConeProperty::DefaultMode) = "DefaultMode";
        CPN.at(QConeProperty::DualMode) = "DualMode";
        CPN.at(QConeProperty::KeepOrder) = "KeepOrder";
        CPN.at(QConeProperty::IntegerHull) = "IntegerHull";
        CPN.at(QConeProperty::ProjectCone) = "ProjectCone";
        CPN.at(QConeProperty::MaximalSubspace) = "MaximalSubspace";
        CPN.at(QConeProperty::ConeDecomposition) = "ConeDecomposition";
        CPN.at(QConeProperty::HSOP) = "HSOP";
        CPN.at(QConeProperty::NoBottomDec) = "NoBottomDec";        
        CPN.at(QConeProperty::PrimalMode) = "PrimalMode";
        CPN.at(QConeProperty::Symmetrize) = "Symmetrize";
        CPN.at(QConeProperty::NoSymmetrization) = "NoSymmetrization";
        CPN.at(QConeProperty::EmbeddingDim) = "EmbeddingDim";
        CPN.at(QConeProperty::Rank) = "Rank";
        CPN.at(QConeProperty::InternalIndex) = "InternalIndex";
        CPN.at(QConeProperty::IsInhomogeneous) = "IsInhomogeneous";
        CPN.at(QConeProperty::UnitGroupIndex) = "UnitGroupIndex";
        CPN.at(QConeProperty::GradingDenom) = "GradingDenom";
        CPN.at(QConeProperty::Equations) = "Equations";
        CPN.at(QConeProperty::Congruences) = "Congruences";
        CPN.at(QConeProperty::ExternalIndex) = "ExternalIndex";
        CPN.at(QConeProperty::HilbertQuasiPolynomial) = "HilbertQuasiPolynomial";
        CPN.at(QConeProperty::IsTriangulationNested) = "IsTriangulationNested";
        CPN.at(QConeProperty::IsTriangulationPartial) = "IsTriangulationPartial";
        CPN.at(QConeProperty::BigInt) = "BigInt";
        CPN.at(QConeProperty::NoSubdivision) = "NoSubdivision";
        CPN.at(QConeProperty::Projection) = "Projection";
        CPN.at(QConeProperty::ProjectionFloat) = "ProjectionFloat";
        CPN.at(QConeProperty::NoProjection) = "NoProjection";
        CPN.at(QConeProperty::NoNestedTri) = "NoNestedTri";
        CPN.at(QConeProperty::Integral) = "Integral";
        CPN.at(QConeProperty::EuclideanIntegral) = "EuclideanIntegral";
        CPN.at(QConeProperty::VirtualMultiplicity) = "VirtualMultiplicity";
        CPN.at(QConeProperty::WeightedEhrhartSeries) = "WeightedEhrhartSeries";
        CPN.at(QConeProperty::WeightedEhrhartQuasiPolynomial) = "WeightedEhrhartQuasiPolynomial";
        CPN.at(QConeProperty::EhrhartSeries) = "EhrhartSeries";
        CPN.at(QConeProperty::EhrhartQuasiPolynomial) = "EhrhartQuasiPolynomial";
        CPN.at(QConeProperty::IsGorenstein) = "IsGorenstein";
        CPN.at(QConeProperty::NoPeriodBound) = "NoPeriodBound";
        CPN.at(QConeProperty::SCIP) = "SCIP";
        CPN.at(QConeProperty::NoLLL) = "NoLLL";
        CPN.at(QConeProperty::NoRelax) = "NoRelax";
        CPN.at(QConeProperty::GeneratorOfInterior) = "GeneratorOfInterior";
        CPN.at(QConeProperty::ExplicitHilbertSeries) = "ExplicitHilbertSeries";
        CPN.at(QConeProperty::NakedDual) = "NakedDual";
        CPN.at(QConeProperty::Descent) = "Descent";
        CPN.at(QConeProperty::NoDescent) = "NoDescent";
        CPN.at(QConeProperty::NoGradingDenom) = "NoGradingDenom";
        CPN.at(QConeProperty::GradingIsPositive) = "GradingIsPositive";
        
        // detect changes in size of Enum, to remember to update CPN!
        static_assert (QConeProperty::EnumSize == 87,
            "ConeProperties Enum size does not fit! Update cone_property.cpp!");
        // assert all fields contain an non-empty string
        for (size_t i=0;  i<QConeProperty::EnumSize; i++) {
            assert(CPN.at(i).size() > 0);
        }
        return CPN;
    }

    const vector<string>& ConePropertyNames() {
        static const vector<string> CPN(initializeCPN());
        return CPN;
    }
}

bool isConeProperty(QConeProperty::Enum& cp, const std::string& s) {
    const vector<string>& CPN = ConePropertyNames();
    for (size_t i=0; i<QConeProperty::EnumSize; i++) {
        if (CPN[i] == s) {
            cp = static_cast<QConeProperty::Enum>(i);
            return true;
        }
    }
    return false;
}

QConeProperty::Enum toConeProperty(const std::string& s) {
    QConeProperty::Enum cp;
    if (isConeProperty(cp, s)) return cp;
    throw BadInputException("Unknown ConeProperty string \"" + s + "\"");
}

const std::string& toString(QConeProperty::Enum cp) {
    return ConePropertyNames()[cp];
}

/* print it in a nice way */
std::ostream& operator<< (std::ostream& out, const ConeProperties& CP){
    for (size_t i=0; i<QConeProperty::EnumSize; i++) {
        if (CP.CPs.test(i)) out << toString(static_cast<QConeProperty::Enum>(i)) << " ";
    }
    return out;
}


} /* end namespace libQnormaliz */

#ifdef NMZ_MIC_OFFLOAD
#pragma offload_attribute (pop)
#endif
