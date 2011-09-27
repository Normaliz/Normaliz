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

#include <list>
#include "cone.h"

namespace libnormaliz {
using namespace std;

template<typename Integer>
Cone<Integer>::Cone(const vector< vector<Integer> >& Input, InputType input_type) {
	initialize();
	if (!Input.empty()) dim = (*(Input.begin())).size();

	switch (input_type) {
		case Type::integral_closure: prepare_input_type_0(Input); break;
		case Type::normalization:    prepare_input_type_1(Input); break;
		case Type::polytope:         prepare_input_type_2(Input); break;
		case Type::rees_algebra:     prepare_input_type_3(Input); break;
		case Type::hyperplanes:
		  prepare_input_type_456(vector<vector<Integer> >(), vector<vector<Integer> >(), Input);
		  break;
		case Type::equations:
		  prepare_input_type_456(vector<vector<Integer> >(), Input, vector<vector<Integer> >());
		  break;
		case Type::congruences:
		  dim--;
		  prepare_input_type_456(Input, vector<vector<Integer> >(), vector<vector<Integer> >());
		  break;
		case Type::lattice_ideal:    prepare_input_type_10(Input); break;
		case Type::grading:
		  errorOutput() << "Grading as only input not supported!" << endl;
		  // no break, go to default
		default:
		  throw BadInputException();
	}
	if(!BC_set) compose_basis_change(Sublattice_Representation<Integer>(dim));
}

template<typename Integer>
Cone<Integer>::Cone(const multimap< InputType , vector< vector<Integer> > >& multi_input_data) {
	initialize();
	
	typename multimap< InputType , vector< vector<Integer> > >::const_iterator it = multi_input_data.begin();
	for(; it != multi_input_data.end(); ++it) {
		if (it->second.size() > 0) {
			dim = it->second.begin()->size();
			if (it->first == Type::congruences) {
				dim--; //congruences have one extra column
			}
			break;
		}
	}
	Matrix<Integer> Inequalities(0,dim), Equations(0,dim), Congruences(0,dim+1);
	for (; it != multi_input_data.end(); ++it) {
		if (it->second.size() == 0) {
			continue;
		}
		switch (it->first) {
			case Type::hyperplanes:
				if (it->second.begin()->size() != dim) {
					errorOutput() << "Dimensions of hyperplanes ("<<it->second.begin()->size()<<") do not match dimension of other constraints ("<<dim<<")!"<<endl;
					throw BadInputException();
				}
				Inequalities.append(it->second);
				break;
			case Type::equations:
				if (it->second.begin()->size() != dim) {
					errorOutput() << "Dimensions of equations ("<<it->second.begin()->size()<<") do not match dimension of other constraints ("<<dim<<")!"<<endl;
					throw BadInputException();
				}
				Equations.append(it->second);
				break;
			case Type::congruences:
				if (it->second.begin()->size() != dim+1) {
					errorOutput() << "Dimensions of congruences ("<<it->second.begin()->size()<<") do not match dimension of other constraints ("<<dim<<")!"<<endl;
					throw BadInputException();
				}
				Congruences.append(it->second);
				break;
			default:
				errorOutput() << "This InputType combination is currently not supported!" << endl;
				throw BadInputException();
		}
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
Sublattice_Representation<Integer> Cone<Integer>::getBasisChange() const{
	return BasisChange;
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getGeneratorsOfToricRing() const {
	return GeneratorsOfToricRing;
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getGenerators() const {
	return Generators;
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getExtremeRays() const {
	return Matrix<Integer>(Generators).submatrix(ExtremeRays).get_elements();
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getSupportHyperplanes() const {
   return SupportHyperplanes;
}

//TODO gehts nicht auch in der SR?
template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getEquations() const {
	size_t rank = BasisChange.get_rank();
	size_t nr_of_equ = dim-rank;
	Matrix<Integer> Equations(nr_of_equ,dim);
	if (nr_of_equ > 0) {
		Lineare_Transformation<Integer> NewLT = Transformation(Matrix<Integer>(getExtremeRays()));
		Matrix<Integer> Help = NewLT.get_right().transpose();
		for (size_t i = 1+rank; i <= dim; i++) {
			Equations.write(i-rank,Help.read(i));
		}
	}
	return Equations.get_elements();
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getCongruences() const {
	return BasisChange.get_congruences().get_elements();
}

template<typename Integer>
multimap< InputType , vector< vector<Integer> > > Cone<Integer>::getConstraints () const {
	multimap<InputType, vector< vector<Integer> > > c;
	c.insert(pair< InputType,vector< vector<Integer> > >(Type::hyperplanes,SupportHyperplanes));
	c.insert(pair< InputType,vector< vector<Integer> > >(Type::equations,getEquations()));
	c.insert(pair< InputType,vector< vector<Integer> > >(Type::congruences,getCongruences()));
	return c;
}


template<typename Integer>
vector< pair<vector<size_t>,Integer> > Cone<Integer>::getTriangulation() const {
	return Triangulation;
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getHilbertBasis() const {
	return HilbertBasis;
}

template<typename Integer>
vector< vector<Integer> > Cone<Integer>::getHt1Elements() const {
	return Ht1Elements;
}

template<typename Integer>
vector<long64> Cone<Integer>::getHVector() const {
	return HSeries.getNumerator();
}

template<typename Integer>
vector<mpz_class> Cone<Integer>::getHilbertPolynomial() {
	if (HSeries.getHilbertQuasiPolynomial().size()==1) {
		return HSeries.getHilbertQuasiPolynomial()[0];
	} else {
		//TODO don't know what to do
		return vector<mpz_class>();
	}
}

template<typename Integer>
vector< vector<mpz_class> > Cone<Integer>::getHilbertQuasiPolynomial() {
	return HSeries.getHilbertQuasiPolynomial();
}

template<typename Integer>
vector<Integer> Cone<Integer>::getLinearForm() const {
	return LinearForm;
}

template<typename Integer>
Integer Cone<Integer>::getMultiplicity() const {
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
Integer Cone<Integer>::getReesPrimaryMultiplicity() const {
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
void Cone<Integer>::prepare_input_type_0(const vector< vector<Integer> >& Input) {
	Generators = Input;
	is_Computed.set(ConeProperty::Generators);

	Sublattice_Representation<Integer> Basis_Change(Matrix<Integer>(Input),true);
	compose_basis_change(Basis_Change);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_1(const vector< vector<Integer> >& Input) {
	Generators = Input;
	is_Computed.set(ConeProperty::Generators);

	Sublattice_Representation<Integer> Basis_Change(Matrix<Integer>(Input),false);
	compose_basis_change(Basis_Change);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_2(const vector< vector<Integer> >& Input) {
	size_t j;
	size_t nr = Input.size();
	if (nr == 0) {
		Generators = Input;
	} else { //append a column of 1
		Generators = vector< vector<Integer> >(nr);
		typename vector< vector<Integer> >::const_iterator it=Input.begin();
		vector<Integer> row(dim+1);
		row[dim]=1;
		for (size_t i=0; i<nr; i++) {
			for (j=0; j<dim; j++) row[j]=Input[i][j];
			Generators[i]=row;
		}
		dim++;
	}
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Matrix<Integer>(Generators),true));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_3(const vector< vector<Integer> >& InputV) {
	Matrix<Integer> Input(InputV);
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
		for(i=1; i<=nr_rows && Prim_Test.read(i,j)==0; i++) {}
		if (i>nr_rows) {
			rees_primary=false;
			break;
		}
	}
	is_Computed.set(ConeProperty::ReesPrimary);
	Generators = Full_Cone_Generators.get_elements();
	dim = Generators.begin()->size();
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Full_Cone_Generators.nr_of_columns()));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_456(const Matrix<Integer>& Congruences, const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities) {

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
void Cone<Integer>::prepare_input_type_45(const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities) {

	// use positive orthant if no inequalities are given
	if (Inequalities.nr_of_rows() == 0) {
		SupportHyperplanes = (Matrix<Integer>(dim)).get_elements();
	} else {
		SupportHyperplanes = Inequalities.get_elements();
	}
	is_Computed.set(ConeProperty::SupportHyperplanes);


	size_t i,j;
	if (Equations.nr_of_rows()>0) {
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
void Cone<Integer>::prepare_input_type_10(const vector< vector<Integer> >& BinomialsV) {
	Matrix<Integer> Binomials(BinomialsV);
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
	GeneratorsOfToricRing = Positive_Embedded_Generators.get_elements();
	is_Computed.set(ConeProperty::GeneratorsOfToricRing);
	dim = Positive_Embedded_Generators.nr_of_columns();
	prepare_input_type_1(GeneratorsOfToricRing);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::setLinearForm (vector<Integer> lf) {
	if (lf.size() != dim) {
		errorOutput() << "Linear form has wrong dimension " << lf.size()
		              << " (should be" << dim << ")" << endl;
		throw BadInputException();
	}
	//check if the linear forms are the same
	if (isComputed(ConeProperty::Generators) && LinearForm == lf) {
		return;
	}
	if (isComputed(ConeProperty::Generators)) {
		vector<Integer> degrees = Matrix<Integer>(Generators).MxV(lf);
		for (size_t i=0; i<degrees.size(); ++i) {
			if (degrees[i]<1) {
				errorOutput() << "Linear form gives non-positive value " << degrees[i]
				              << " for generator " << i+1 << "." << endl;
				throw BadInputException();
			}
		}
	}
	LinearForm = lf;
	is_Computed.set(ConeProperty::LinearForm);

	//remove data that depends on the grading 
	Ht1Elements.clear();
    HilbertQuasiPolynomial.clear();
	is_Computed.reset(ConeProperty::IsHt1Generated);
	is_Computed.reset(ConeProperty::IsHt1ExtremeRays);
	is_Computed.reset(ConeProperty::IsHt1HilbertBasis);
	is_Computed.reset(ConeProperty::Ht1Elements);
	is_Computed.reset(ConeProperty::HVector);
	is_Computed.reset(ConeProperty::HilbertPolynomial);
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
			compute(Mode::hilbertBasisPolynomial);
		} else {
			compute(Mode::hilbertPolynomial);
		}
	} else { //no H-Vector
		if(ToCompute.test(ConeProperty::HilbertBasis)) {
			if(ToCompute.test(ConeProperty::Triangulation)) {
				compute(Mode::hilbertBasisTriangulation);
			} else {
				compute(Mode::hilbertBasisLarge);
			}
		} else { //no Hilbert basis
			if(ToCompute.test(ConeProperty::Triangulation)) {
				compute(Mode::volumeTriangulation);
				if(ToCompute.test(ConeProperty::Ht1Elements)) {
					compute(Mode::height1Elements);
				}
			} else { //no triangulation
				if(ToCompute.test(ConeProperty::Ht1Elements)) {
					compute(Mode::height1Elements);
				} else if(ToCompute.test(ConeProperty::SupportHyperplanes)) {
					compute(Mode::supportHyperplanes);
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
void Cone<Integer>::compute(ComputationMode mode) {
	if (mode == Mode::dual) {
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
			Generators = BasisChange.from_sublattice(Extreme_Rays).get_elements();
			//sort Generators to get deterministic triangulations
			sort (Generators.begin(), Generators.end());
			is_Computed.set(ConeProperty::Generators);
			//get minmal set of support_hyperplanes
			Matrix<Integer> Supp_Hyp = Dual_Cone.getGenerators().submatrix(Dual_Cone.getExtremeRays());
			SupportHyperplanes = BasisChange.from_sublattice_dual(Supp_Hyp).get_elements();

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

	// Create a Full_Cone FC
	Full_Cone<Integer> FC(BasisChange.to_sublattice(Matrix<Integer>(Generators)));

	// Give extra data to FC
	if ( isComputed(ConeProperty::LinearForm) ) {
		FC.Linear_Form = BasisChange.to_sublattice_dual(LinearForm);
		FC.is_Computed.set(ConeProperty::LinearForm);
	}

	// Start computations in FC
	switch (mode) {
	case Mode::hilbertBasisTriangulation:
		FC.triangulation_hilbert_basis();
		break;
	case Mode::hilbertBasisLarge:
		FC.hilbert_basis();
		break;
	case Mode::supportHyperplanes:
		// workaround for not dualizing twice
		if (isComputed(ConeProperty::Generators)
		 && isComputed(ConeProperty::SupportHyperplanes)) {
			vector< vector<Integer> > vvSH = BasisChange.to_sublattice_dual(Matrix<Integer>(SupportHyperplanes)).get_elements();
			FC.Support_Hyperplanes = list< vector<Integer> >(vvSH.begin(), vvSH.end());
			FC.is_Computed.set(ConeProperty::SupportHyperplanes);
		}
		FC.support_hyperplanes();
		break;
	case Mode::volumeTriangulation:
		FC.support_hyperplanes_triangulation();
		break;
	case Mode::volumeLarge:
		FC.support_hyperplanes_triangulation_pyramid();
		break;
	case Mode::height1Elements:
		FC.ht1_elements();
		break;
	case Mode::hilbertPolynomial:
		FC.hilbert_polynomial();
		break;
	case Mode::hilbertPolynomialLarge:
		FC.hilbert_polynomial_pyramid();
		break;
	case Mode::hilbertBasisPolynomial:
		FC.hilbert_basis_polynomial();
		break;
	case Mode::hilbertBasisPolynomialLarge:
		FC.hilbert_basis_polynomial_pyramid();
		break;
	default: //should not happen
		errorOutput()<<"Unknown computation mode: \""<<static_cast<int>(mode)<<"\"!"<<endl;
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
	
	if (rees_primary && FC.isComputed(ConeProperty::Triangulation)) {
		//here are some computations involved, made first so that data can be deleted in FC later
		ReesPrimaryMultiplicity = FC.primary_multiplicity();
		is_Computed.set(ConeProperty::ReesPrimaryMultiplicity);
	}
	
	if (FC.isComputed(ConeProperty::Generators)) {
		Generators = BasisChange.from_sublattice(FC.getGenerators()).get_elements();
		is_Computed.set(ConeProperty::Generators);
	}
	if (FC.isComputed(ConeProperty::ExtremeRays)) {
		ExtremeRays = FC.getExtremeRays();
		is_Computed.set(ConeProperty::ExtremeRays);
	}
	if (FC.isComputed(ConeProperty::SupportHyperplanes)) {
		SupportHyperplanes = BasisChange.from_sublattice_dual(FC.getSupportHyperplanes()).get_elements();
		is_Computed.set(ConeProperty::SupportHyperplanes);
	}
	if (FC.isComputed(ConeProperty::Triangulation)) {
		size_t tri_size = FC.Triangulation.size();
		Triangulation = vector< pair<vector<size_t>, Integer> >();
		Triangulation.reserve(tri_size);
		for (size_t i = 0; i<tri_size; ++i) {
			Triangulation.push_back(FC.Triangulation.front());
			FC.Triangulation.pop_front();
		}
		is_Computed.set(ConeProperty::Triangulation);
	}
	if (FC.isComputed(ConeProperty::Multiplicity)) {
		multiplicity = FC.getMultiplicity();
		is_Computed.set(ConeProperty::Multiplicity);
	}
	if (FC.isComputed(ConeProperty::HilbertBasis)) {
		HilbertBasis = BasisChange.from_sublattice(FC.getHilbertBasis()).get_elements();
		is_Computed.set(ConeProperty::HilbertBasis);
	}
	if (FC.isComputed(ConeProperty::Ht1Elements)) {
		Ht1Elements = BasisChange.from_sublattice(FC.getHt1Elements()).get_elements();
		is_Computed.set(ConeProperty::Ht1Elements);
	}
	if (FC.isComputed(ConeProperty::HVector)) {
		//TODO HilbertSeries
		HSeries = FC.Hilbert_Series;
		is_Computed.set(ConeProperty::HVector);
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

	if (verbose) {
		verboseOutput() << " done." <<endl;
	}
}

} // end namespace libnormaliz
