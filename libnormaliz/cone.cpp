/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#include <map>
#include "cone.h"

namespace libnormaliz {
using namespace std;

template<typename T> void list_vector_print(const list< vector<T> >& l) {
	typename list< vector<T> >::const_iterator it = l.begin();
	for (; it != l.end(); ++it) {
		v_read(*it);
	}
	verboseOutput()<<endl;
}

template<typename Integer>
Cone<Integer>::Cone(const list< vector<Integer> >& Input, int input_type) {
	initialize();
	if (!Input.empty()) dim = (*(Input.begin())).size();

	switch (input_type){
		case  0: prepare_input_type_0(Input); break;
		case  1: prepare_input_type_1(Input); break;
		case  2: prepare_input_type_2(Input); break;
		case  3: prepare_input_type_3(Input); break;
		case  4: prepare_input_type_456(list<vector<Integer> >(), list<vector<Integer> >(), Input); break;
		case  5: prepare_input_type_456(list<vector<Integer> >(), Input, list<vector<Integer> >()); break;
		case  6: dim--; prepare_input_type_456(Input, list<vector<Integer> >(), list<vector<Integer> >()); break;
		case 10: prepare_input_type_10(Input); break;
		default: throw input_type; //TODO make a good exception
	}
	if(!BC_set) compose_basis_change(Sublattice_Representation<Integer>(dim));
}

template<typename Integer>
Cone<Integer>::Cone(const list< vector<Integer> >& Inequalities, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Congruences) {
	initialize();
	list<size_t> dimensions = list<size_t>();
	if (Inequalities.size()>0) dimensions.push_back(Inequalities.begin()->size());
	if (Equations.size()>0) dimensions.push_back(Equations.begin()->size());
	if (Congruences.size()>0) dimensions.push_back(Congruences.begin()->size()-1);
	if (!dimensions.empty()) {
		dim = *dimensions.begin();
		dimensions.pop_front();
	}
	while(!dimensions.empty()) {
		if (dim != *dimensions.begin()) {
			errorOutput() << "Error: dimensions of input matrices do not match!";
			throw NormalizException();
		}
		dimensions.pop_front();
	}
	if(!BC_set) compose_basis_change(Sublattice_Representation<Integer>(dim));
	prepare_input_type_456(Congruences, Equations, Inequalities);
}

/* only used by the constructors */
template<typename Integer>
void Cone<Integer>::initialize() {
	BC_set=false;
	is_Computed = bitset<ConeProperty::EnumSize>();  //initialized to false
	dim = 0;
	rees_primary = false;
}


/* check what is computed */
template<typename Integer>
bool Cone<Integer>::isComputed(ConeProperty::Enum prop) const {
	return is_Computed.test(prop);
}


/* getter */
template<typename Integer>
Sublattice_Representation<Integer> const& Cone<Integer>::getBasisChange() const{
	return BasisChange;
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getGeneratorsOfToricRing() const {
	return GeneratorsOfToricRing;
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getGenerators() const {
	return Generators;
}

template<typename Integer>
vector<bool> const& Cone<Integer>::getExtremeRays() const {
	return ExtremeRays;
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getSupportHyperplanes() const {
   return SupportHyperplanes;
}

template<typename Integer>
list< pair<vector<size_t>,Integer> > const& Cone<Integer>::getTriangulation() const {
	return Triangulation;
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHilbertBasis() const {
	return HilbertBasis;
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHt1Elements() const {
	return Ht1Elements;
}

template<typename Integer>
vector<Integer> const& Cone<Integer>::getHVector() const {
	return HVector;
}

template<typename Integer>
vector<Integer> const& Cone<Integer>::getHilbertPolynomial() const {
	return HilbertPolynomial;
}

template<typename Integer>
vector<Integer> const& Cone<Integer>::getLinearForm() const {
	return LinearForm;
}

template<typename Integer>
Integer const& Cone<Integer>::getMultiplicity() const {
	return multiplicity;
}

template<typename Integer>
bool Cone<Integer>::isPointed() const {
	return pointed;
}

template<typename Integer>
bool Cone<Integer>::isHt1ExtremeRays() const {
	return ht1_extreme_rays;
}

template<typename Integer>
bool Cone<Integer>::isHt1HilbertBasis() const {
	return ht1_hilbert_basis;
}

template<typename Integer>
bool Cone<Integer>::isIntegrallyClosed() const {
	return integrally_closed;
}

template<typename Integer>
bool Cone<Integer>::isReesPrimary() const {
	return rees_primary;
}

template<typename Integer>
Integer const& Cone<Integer>::getReesPrimaryMultiplicity() const {
	return ReesPrimaryMultiplicity;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::compose_basis_change(const Sublattice_Representation<Integer>& BC) {
	if (BC_set) {
		BasisChange.compose(BC);
	} else {
		BasisChange = BC;
		BC_set = true;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_0(const list< vector<Integer> >& Input) {
	Generators = Input;
	is_Computed.set(ConeProperty::Generators);

	Sublattice_Representation<Integer> Basis_Change(Matrix<Integer>(Input),true);
	compose_basis_change(Basis_Change);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_1(const list< vector<Integer> >& Input) {
	Generators = Input;
	is_Computed.set(ConeProperty::Generators);

	Sublattice_Representation<Integer> Basis_Change(Matrix<Integer>(Input),false);
	compose_basis_change(Basis_Change);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_2(const list< vector<Integer> >& Input) {
	size_t j;
	size_t nr = Input.size();
	if (nr == 0) {
		Generators = Input;
	} else { //append a column of 1
		Generators = list< vector<Integer> >();
		typename list< vector<Integer> >::const_iterator it=Input.begin();
		vector<Integer> row(dim+1);
		row[dim]=1;
		for (; it!=Input.end(); ++it) {
			for (j=0; j<dim; j++) row[j]=(*it)[j];
			Generators.push_back(row);
		}
		dim++;
	}
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Matrix<Integer>(Generators),true));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_3(const list< vector<Integer> >& InputL) {
	Matrix<Integer> Input(InputL);  //TODO handle it better
	int i,j,k,l,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	rees_primary=true;
	Integer number;
	Matrix<Integer> Full_Cone_Generators(nr_rows+nr_columns,nr_columns+1,0);
	for (i = 1; i <= nr_columns; i++) {
		Full_Cone_Generators.write(i,i,1);
	}
	for(i=1; i<=nr_rows; i++){
		Full_Cone_Generators.write(i+nr_columns,nr_columns+1,1);
		for(j=1; j<=nr_columns; j++) {
			number=Input.read(i,j);
			Full_Cone_Generators.write(i+nr_columns,j,number);
		}
	}
	Matrix<Integer> Prim_Test=Input;
	for(i=1; i<=nr_rows; i++){           //preparing the  matrix for primarity test
		k=0;
		for(j=1; j<=nr_columns; j++) {
			if (k<2) {
				if (Input.read(i,j)!=0 )
					k++;
			}
			if (k==2) {
				for (l = 1; l <= nr_columns; l++) {
					Prim_Test.write(i,l,0);
				}
				break;
			}
		}
	}
	for(j=1; j<=nr_columns; j++){         //primarity test
		for(i=1; i<=nr_rows && Prim_Test.read(i,j)==0; i++);
		if (i>nr_rows) {
			rees_primary=false;
			is_Computed.set(ConeProperty::ReesPrimary);
			break;
		}
	}
	Generators = Full_Cone_Generators.to_list();
	dim = Generators.begin()->size();
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Full_Cone_Generators.nr_of_columns()));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_456(const list< vector<Integer> >& CongruencesL, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Inequalities) {
	Matrix<Integer> Congruences(CongruencesL); //TODO handle it better

	size_t nr_cong = Congruences.nr_of_rows();
	// handle Congurences
	if (nr_cong > 0) {
		size_t i,j;

		//add slack variables
		Matrix<Integer> Cong_Slack(nr_cong, dim+nr_cong);
		for (i = 1; i <= nr_cong; i++) {
			for (j = 1; j <= dim; j++) {
				Cong_Slack.write(i,j,Congruences.read(i,j));
			}
			Cong_Slack.write(i,dim+i,Congruences.read(i,dim+1));
		}

		//compute kernel
		Lineare_Transformation<Integer> Diagonalization = Transformation(Cong_Slack);
		size_t rank = Diagonalization.get_rank();
		Matrix<Integer> H = Diagonalization.get_right();
		Matrix<Integer> Ker_Basis_Transpose(dim, dim+nr_cong-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim+nr_cong; j++) {
				Ker_Basis_Transpose.write(i,j-rank,H.read(i,j));
			}
		}

		//TODO now a new linear transformation is computed, necessary??
		Sublattice_Representation<Integer> Basis_Change(Ker_Basis_Transpose.transpose(),false);
		compose_basis_change(Basis_Change);
	}

	prepare_input_type_45(Equations, Inequalities);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_45(const list< vector<Integer> >& Equations, const list< vector<Integer> >& Inequalities) {

	// use positive orthant if no inequalities are given
	if (Inequalities.size() == 0) {
		SupportHyperplanes = (Matrix<Integer>(dim)).to_list();
	} else {
		SupportHyperplanes = Inequalities;
	}
	is_Computed.set(ConeProperty::SupportHyperplanes);


	size_t i,j;
	if (Equations.size()>0) {
		Lineare_Transformation<Integer> Diagonalization = Transformation(BasisChange.to_sublattice_dual(Equations));
		size_t rank=Diagonalization.get_rank();

		Matrix<Integer> Help=Diagonalization.get_right();
		Matrix<Integer> Ker_Basis_Transpose(dim,dim-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim; j++) {
				Ker_Basis_Transpose.write(i,j-rank,Help.read(i,j));
			}
		}
		Sublattice_Representation<Integer> Basis_Change(Ker_Basis_Transpose.transpose(),true);
		compose_basis_change(Basis_Change);
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_10(const list< vector<Integer> >& BinomialsL) {
	Matrix<Integer> Binomials(BinomialsL); //TODO geschickter machen
	size_t i,j, nr_of_monoid_generators = dim;
	Lineare_Transformation<Integer> Diagonalization=Transformation(Binomials);
	size_t rank=Diagonalization.get_rank();
	Matrix<Integer> Help=Diagonalization.get_right();
	Matrix<Integer> Generators(nr_of_monoid_generators,nr_of_monoid_generators-rank);
	for (i = 1; i <= nr_of_monoid_generators; i++) {
		for (j = rank+1; j <= nr_of_monoid_generators; j++) {
			Generators.write(i,j-rank,Help.read(i,j));
		}
	}
	Full_Cone<Integer> FC(Generators);
	//TODO verboseOutput(), what is happening here?
	FC.support_hyperplanes();
	Matrix<Integer> Supp_Hyp=FC.getSupportHyperplanes();
	Matrix<Integer> Selected_Supp_Hyp_Trans=(Supp_Hyp.submatrix(Supp_Hyp.max_rank_submatrix_lex())).transpose();
	Matrix<Integer> Positive_Embedded_Generators=Generators.multiplication(Selected_Supp_Hyp_Trans);
	GeneratorsOfToricRing = Positive_Embedded_Generators.to_list();
	is_Computed.set(ConeProperty::GeneratorsOfToricRing);
	dim = Positive_Embedded_Generators.nr_of_columns();
	prepare_input_type_1(GeneratorsOfToricRing);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::compute(ConeProperties ToCompute) {
	ToCompute.reset(is_Computed); // already computed

	/* add preconditions */
	if(ToCompute.test(ConeProperty::Multiplicity))       ToCompute.set(ConeProperty::Triangulation);
	if(ToCompute.test(ConeProperty::IsIntegrallyClosed)) ToCompute.set(ConeProperty::HilbertBasis);
	if(ToCompute.test(ConeProperty::IsHt1HilbertBasis))  ToCompute.set(ConeProperty::HilbertBasis);
	if(ToCompute.test(ConeProperty::IsHt1ExtremeRays))   ToCompute.set(ConeProperty::ExtremeRays);
	if(ToCompute.test(ConeProperty::LinearForm))         ToCompute.set(ConeProperty::ExtremeRays);
	if(ToCompute.test(ConeProperty::ExtremeRays))        ToCompute.set(ConeProperty::SupportHyperplanes);
	if(ToCompute.test(ConeProperty::IsPointed))          ToCompute.set(ConeProperty::SupportHyperplanes);
	if(ToCompute.test(ConeProperty::HilbertPolynomial))  ToCompute.set(ConeProperty::HVector);


	/* find correct mode */
	if (ToCompute.test(ConeProperty::HVector) ) {
		if(ToCompute.test(ConeProperty::HilbertBasis)) {
			compute("hilbert_basis_polynomial");
		} else {
			compute("hilbert_polynomial");
		}
	} else { //no H-Vector
		if(ToCompute.test(ConeProperty::HilbertBasis)) {
			if(ToCompute.test(ConeProperty::Triangulation)) {
				compute("triangulation_hilbert_basis");
			} else {
				compute("hilbert_basis");
			}
		} else { //no Hilbert basis
			if(ToCompute.test(ConeProperty::Triangulation)) {
				compute("triangulation");
				if(ToCompute.test(ConeProperty::Ht1Elements)) {
					compute("ht1_elements");
				}
			} else { //no triangulation
				if(ToCompute.test(ConeProperty::Ht1Elements)) {
					compute("ht1_elements");
				} else if(ToCompute.test(ConeProperty::SupportHyperplanes)) {
					compute("support_hyperplanes");
				}
			}
		}
	}

	/* check if everything is computed*/
	ToCompute.reset(is_Computed); //remove what is now computed
	if (ToCompute.any()) {
		errorOutput() << "Warning: Cone could not compute everything, that it was asked for!"<<endl;
		errorOutput() << "Missing: "; ToCompute.print(errorOutput());
	}
}



template<typename Integer>
void Cone<Integer>::compute(const string& computation_type) {
	if (computation_type == "dual") {
		compute_dual();
		return;
	}

	//create Generators from SupportHyperplanes
	if (!isComputed(ConeProperty::Generators) && isComputed(ConeProperty::SupportHyperplanes)) {
		if (verbose) {
			verboseOutput() <<endl<< "Computing extreme rays as support hyperplanes of the dual cone:";
		}
		Full_Cone<Integer> Dual_Cone(BasisChange.to_sublattice_dual(Matrix<Integer>(SupportHyperplanes)));
		Dual_Cone.support_hyperplanes();
		if (Dual_Cone.isComputed(ConeProperty::SupportHyperplanes)) {
			//get the extreme rays of the primal cone
			Matrix<Integer> Extreme_Rays=Dual_Cone.getSupportHyperplanes();
			Generators = BasisChange.from_sublattice(Extreme_Rays).to_list();
			//sort Generators to get deterministic triangulations
			Generators.sort();
			is_Computed.set(ConeProperty::Generators);
			//get minmal set of support_hyperplanes
			Matrix<Integer> Supp_Hyp = Dual_Cone.getGenerators().submatrix(Dual_Cone.getExtremeRays());
			SupportHyperplanes = BasisChange.from_sublattice_dual(Supp_Hyp).to_list();

			Sublattice_Representation<Integer> Basis_Change(Extreme_Rays,true);
			compose_basis_change(Basis_Change);
		}
	}

	if (!isComputed(ConeProperty::Generators)) {
		errorOutput()<<"FATAL ERROR: Could not get Generators. This should not happen!"<<endl;
		throw NormalizException();
	}

	//TODO workaround for zero cones :((
	//die dimension bleibt in der liste nicht gespeichert, wenn sie leer ist, darum passt es dann beim transformieren nicht
	if(Generators.size()==0) {
		return;
	}

	Full_Cone<Integer> FC(BasisChange.to_sublattice(Matrix<Integer>(Generators)));

	if (computation_type=="triangulation_hilbert_basis") {
		FC.triangulation_hilbert_basis();
	} else if (computation_type=="hilbert_basis") {
		FC.hilbert_basis();
	} else if (computation_type=="support_hyperplanes" ||
	           computation_type=="support_hyperplanes_pyramid") {
		if (isComputed(ConeProperty::Generators) && isComputed(ConeProperty::SupportHyperplanes)) {
			//this is the workaround for not dualizing twice
			FC.Support_Hyperplanes=BasisChange.to_sublattice_dual(Matrix<Integer>(SupportHyperplanes)).to_list();
			FC.is_Computed.set(ConeProperty::SupportHyperplanes);
		}
		FC.support_hyperplanes();
	} else if (computation_type=="triangulation") {
		FC.support_hyperplanes_triangulation();
	} else if (computation_type=="volume") {
		FC.support_hyperplanes_triangulation_pyramid();
	} else if (computation_type=="ht1_elements") {
		FC.ht1_elements();
	} else if (computation_type=="hilbert_polynomial") {
		FC.hilbert_polynomial();
	} else if (computation_type=="hilbert_polynomial_pyramid") {
		FC.hilbert_polynomial_pyramid();
	} else if (computation_type=="hilbert_basis_polynomial") {
		FC.hilbert_basis_polynomial();
	} else if (computation_type=="hilbert_basis_polynomial_pyramid") {
		FC.hilbert_basis_polynomial_pyramid();
	} else {
		errorOutput()<<"Unknown computation_type: \""<<computation_type<<"\"!"<<endl;
		throw NormalizException();
	}
	extract_data(FC);
}

template<typename Integer>
void Cone<Integer>::compute_dual() {
	if(isComputed(ConeProperty::Generators) && !isComputed(ConeProperty::SupportHyperplanes)){
		if (verbose) {
			verboseOutput() <<endl<< "Computing support hyperplanes for the dual mode:";
		}
		Full_Cone<Integer> Tmp_Cone(BasisChange.to_sublattice(Matrix<Integer>(Generators)));
		Tmp_Cone.support_hyperplanes();
		extract_data(Tmp_Cone);
	}

	if (!isComputed(ConeProperty::SupportHyperplanes)) {
		errorOutput()<<"FATAL ERROR: Could not get SupportHyperplanes. This should not happen!"<<endl;
		throw NormalizException();
	}

	size_t i,j;
	Matrix<Integer> Inequ_on_Ker = BasisChange.to_sublattice_dual(Matrix<Integer>(SupportHyperplanes));
	size_t newdim = Inequ_on_Ker.nr_of_columns();
	Integer norm;
	vector< Integer > hyperplane;
	multimap <Integer , vector <Integer> >  Help;
	typename multimap <Integer , vector <Integer> >::const_iterator ii;
	for (i = 1; i <= Inequ_on_Ker.nr_of_rows() ; i++) {
		hyperplane=Inequ_on_Ker.read(i);
		norm=0;
		for (j = 0; j <newdim; j++) {
			norm+=Iabs(hyperplane[j]);
		}
		Help.insert(pair <Integer , vector <Integer> > (norm,hyperplane));
	}
	Matrix<Integer> Equations_Ordered(Inequ_on_Ker.nr_of_rows(),newdim);
	i=1;
	for (ii=Help.begin(); ii != Help.end(); ii++) {
		Equations_Ordered.write(i,(*ii).second);
		i++;
	}
	Cone_Dual_Mode<Integer> ConeDM(Equations_Ordered);
	ConeDM.hilbert_basis_dual();
	//ConeDM zu einem Full_Cone<Integer> machen
	if ( ConeDM.Generators.rank() < ConeDM.dim ) {
		Sublattice_Representation<Integer> SR(ConeDM.Generators,true);
		ConeDM.to_sublattice(SR);
		compose_basis_change(SR);
	}
	Full_Cone<Integer> FC(ConeDM);
	FC.dual_mode();
	extract_data(FC);
}

template<typename Integer>
void Cone<Integer>::extract_data(Full_Cone<Integer>& FC) {
	//this function extracts ALL available data from the Full_Cone
	//even if it was in Cone already <- this may change
	//it is possible to delete the data in Full_Cone after extracting it

	if(verbose) {
		verboseOutput() << "transforming data...";
	}
	
	if (rees_primary) {
		//here are some computations involved, made first so that data can be deleted in FC later
		ReesPrimaryMultiplicity = FC.primary_multiplicity();
		is_Computed.set(ConeProperty::ReesPrimary);
	}
	
	if (FC.isComputed(ConeProperty::Generators)) {
		Generators = BasisChange.from_sublattice(FC.getGenerators()).to_list();
		is_Computed.set(ConeProperty::Generators);
	}
	if (FC.isComputed(ConeProperty::ExtremeRays)) {
		ExtremeRays = FC.getExtremeRays();
		is_Computed.set(ConeProperty::ExtremeRays);
	}
	if (FC.isComputed(ConeProperty::SupportHyperplanes)) {
		SupportHyperplanes = BasisChange.from_sublattice_dual(FC.getSupportHyperplanes()).to_list();
		is_Computed.set(ConeProperty::SupportHyperplanes);
	}
	if (FC.isComputed(ConeProperty::Triangulation)) {
		Triangulation.clear();
		Triangulation.splice(Triangulation.begin(),FC.Triangulation);
		is_Computed.set(ConeProperty::Triangulation);
	}
	if (FC.isComputed(ConeProperty::Multiplicity)) {
		multiplicity = FC.getMultiplicity();
		is_Computed.set(ConeProperty::Multiplicity);
	}
	if (FC.isComputed(ConeProperty::HilbertBasis)) {
		HilbertBasis = BasisChange.from_sublattice(FC.getHilbertBasis()).to_list();
		is_Computed.set(ConeProperty::HilbertBasis);
	}
	if (FC.isComputed(ConeProperty::Ht1Elements)) {
		Ht1Elements = BasisChange.from_sublattice(FC.getHt1Elements()).to_list();
		is_Computed.set(ConeProperty::Ht1Elements);
	}
	if (FC.isComputed(ConeProperty::HVector)) {
		HVector = FC.getHVector();
		is_Computed.set(ConeProperty::HVector);
	}
	if (FC.isComputed(ConeProperty::HilbertPolynomial)) {
		HilbertPolynomial = FC.getHilbertPolynomial();
		is_Computed.set(ConeProperty::HilbertPolynomial);
	}
	if (FC.isComputed(ConeProperty::IsPointed)) {
		pointed = FC.isPointed();
		is_Computed.set(ConeProperty::IsPointed);
	}
	if (FC.isComputed(ConeProperty::IsHt1ExtremeRays)) {
		ht1_extreme_rays = FC.isHt1ExtremeRays();
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
	}
	if (FC.isComputed(ConeProperty::LinearForm)) {
		LinearForm = BasisChange.from_sublattice_dual(FC.getLinearForm());
		is_Computed.set(ConeProperty::LinearForm);
	}
	if (FC.isComputed(ConeProperty::IsHt1HilbertBasis)) {
		ht1_hilbert_basis = FC.isHt1HilbertBasis();
		is_Computed.set(ConeProperty::IsHt1HilbertBasis);
	}
	if (FC.isComputed(ConeProperty::IsIntegrallyClosed)) {
		integrally_closed = FC.isIntegrallyClosed();
		is_Computed.set(ConeProperty::IsIntegrallyClosed);
	}

	if(verbose) {
		verboseOutput() << " done" <<endl;
	}
}

} // end namespace libnormaliz
