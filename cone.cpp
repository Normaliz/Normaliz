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

#include "cone.h"

namespace libnormaliz {

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
}

template<typename Integer>
Cone<Integer>::Cone(const list< vector<Integer> >& Inequalities, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Congruences) {
	initialize();
	prepare_input_type_456(Congruences, Equations, Inequalities);
}

/* only used by the constructors */
template<typename Integer>
void Cone<Integer>::initialize() {
	BC_set=false; OrigGens_set=false;
	is_Computed =  bitset<ConeProperty::EnumSize>();  //initialized to false
	dim = 0;
}


/* check what is computed */
template<typename Integer>
bool Cone<Integer>::isComputed(ConeProperty::Enum prop) const {
	return is_Computed.test(prop);
}


/* getter */
template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getExtremeRays() const {
	return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getSupportHyperplanes() const {
   return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getTriangulation() const {
	return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHilbertBasis() const {
	return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHt1Elements() const {
	return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHVector() const {
	return Generators; //TODO implement
}

template<typename Integer>
list< vector<Integer> > const& Cone<Integer>::getHilbertPolynomial() const {
	return Generators; //TODO implement
}

template<typename Integer>
vector<Integer> const& Cone<Integer>::getLinearFunction() const {
	return *(Generators.begin()); //TODO implement
}


template<typename Integer>
bool Cone<Integer>::isPointed() const {
	return false; //TODO implement
}

template<typename Integer>
bool Cone<Integer>::isHt1Generated() const {
	return false; //TODO implement
}

template<typename Integer>
bool Cone<Integer>::isHt1ExtremeRays() const {
	return false; //TODO implement
}

template<typename Integer>
bool Cone<Integer>::isHt1HilbertBasis() const {
	return false; //TODO implement
}

template<typename Integer>
bool Cone<Integer>::isIntegrallyClosed() const {
	return false; //TODO implement
}


//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::compose_basis_change(const Sublattice_Representation<Integer>& BC) {
	if (BC_set) {
		ChangeToFullDim.compose(BC);
	} else {
		ChangeToFullDim = BC;
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
	int j;
	int nr = Input.size();
	if (nr == 0) {
		Generators = Input;
	} else { //append a column of 1
		Generators = list< vector<Integer> >();
		typename list< vector<Integer> >::const_iterator it=Input.begin();
		int nc = (*it).size();
		vector<Integer> row(nc+1);
		row[nc]=1;
		for (; it!=Input.end(); ++it) {
			for (j=0; j<nc; j++) row[j]=(*it)[j];
			Generators.push_back(row);
		}
	}
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Matrix<Integer>(Generators),true));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_3(const list< vector<Integer> >& InputL) {
	Matrix<Integer> Input(InputL);  //TODO handle it better
	int i,j,k,l,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	bool primary=true;
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
			primary=false;
			break;
		}
	}
	Generators = Full_Cone_Generators.to_list();
	is_Computed.set(ConeProperty::Generators);

	compose_basis_change(Sublattice_Representation<Integer>(Full_Cone_Generators.nr_of_columns()));
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_456(const list< vector<Integer> >& CongruencesL, const list< vector<Integer> >& EquationsL, const list< vector<Integer> >& InequalitiesL) {
	Matrix<Integer> Congruences(CongruencesL); //TODO handle it better
	Matrix<Integer> Equations(EquationsL);
	Matrix<Integer> Inequalities(InequalitiesL);

	int nr_cong = Congruences.nr_of_rows();
	if (nr_cong > 0) {
		dim = Congruences.nr_of_columns() -1;
		if (Equations.nr_of_rows() > 0 &&  Equations.nr_of_columns() != dim) {
			cerr << "Error: dimensions of input matrices do not match!";
			throw 1; //TODO exception
		}
	} else if (Equations.nr_of_rows() > 0) {
		dim = Equations.nr_of_columns();
	} else if (Inequalities.nr_of_rows() > 0) {
		dim = Inequalities.nr_of_columns();
	}


	// use positive orthant if no inequalities are given
	if (Inequalities.nr_of_rows() == 0) {
		Inequalities = Matrix<Integer>(Equations.nr_of_columns());
	}

	// handle Congurences
	if (nr_cong > 0) {
		int i,j;

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
		int rank = Diagonalization.get_rank();
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
		Equations = Basis_Change.to_sublattice_dual(Equations);
		Inequalities = Basis_Change.to_sublattice_dual(Inequalities);
	}

	prepare_input_type_45(Equations, Inequalities);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_45(const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities) {
	int i,j,dim=Equations.nr_of_columns();
	Lineare_Transformation<Integer> Diagonalization=Transformation(Equations);
	int rank=Diagonalization.get_rank();

	Matrix<Integer> Help=Diagonalization.get_right();
	Matrix<Integer> Ker_Basis_Transpose(dim,dim-rank);
	for (i = 1; i <= dim; i++) {
		for (j = rank+1; j <= dim; j++) {
			Ker_Basis_Transpose.write(i,j-rank,Help.read(i,j));
		}
	}
	Sublattice_Representation<Integer> Basis_Change(Ker_Basis_Transpose.transpose(),true);
	compose_basis_change(Basis_Change);
	Matrix<Integer> Inequ_on_Ker = Basis_Change.to_sublattice_dual(Inequalities);

	//TODO set SH and isComp(SH)
	cerr<<"WAHH! not implemeted yet!"<<endl;
	throw 23;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::prepare_input_type_10(const list< vector<Integer> >& BinomialsL) {
	Matrix<Integer> Binomials(BinomialsL); //TODO geschickter machen
	int i,j, nr_of_monoid_generators = dim;
	Lineare_Transformation<Integer> Diagonalization=Transformation(Binomials);
	int rank=Diagonalization.get_rank();
	Matrix<Integer> Help=Diagonalization.get_right();
	Matrix<Integer> Generators(nr_of_monoid_generators,nr_of_monoid_generators-rank);
	for (i = 1; i <= nr_of_monoid_generators; i++) {
		for (j = rank+1; j <= nr_of_monoid_generators; j++) {
			Generators.write(i,j-rank,Help.read(i,j));
		}
	}
	Full_Cone<Integer> FC(Generators);
	//TODO cout, what is happening here?
	FC.support_hyperplanes();
	Matrix<Integer> Supp_Hyp=FC.read_support_hyperplanes();
	Matrix<Integer> Selected_Supp_Hyp_Trans=(Supp_Hyp.submatrix(Supp_Hyp.max_rank_submatrix_lex())).transpose();
	Matrix<Integer> Positive_Embedded_Generators=Generators.multiplication(Selected_Supp_Hyp_Trans);
	OriginalGenerators = Positive_Embedded_Generators.to_list();
	OrigGens_set=true;
	dim = Positive_Embedded_Generators.nr_of_columns();
	prepare_input_type_1(OriginalGenerators);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone<Integer>::compute(const string& computation_type) {
	//TODO check what is set, the code below should transfer, if necessary
/*	if(computation_type!="dual"){
		if (verbose) {
			cout <<endl<< "Computing extreme rays as support hyperplanes of the dual cone:";
		}
		Full_Cone<Integer> Dual_Cone(Inequ_on_Ker);
		Dual_Cone.support_hyperplanes();
		Matrix<Integer> Extreme_Rays=Dual_Cone.read_support_hyperplanes();
		prepare_input_type_0(Extreme_Rays);
	}
	if(computation_type=="dual"){
		dim = Inequ_on_Ker.nr_of_columns();
		Integer norm;
		vector< Integer > hyperplane;
		multimap <Integer , vector <Integer> >  Help;
		typename multimap <Integer , vector <Integer> >::const_iterator ii;
		for (i = 1; i <= Inequ_on_Ker.nr_of_rows() ; i++) {
			hyperplane=Inequ_on_Ker.read(i);
			norm=0;
			for (j = 0; j <dim; j++) {
				norm+=Iabs(hyperplane[j]);
			}
			Help.insert(pair <Integer , vector <Integer> > (norm,hyperplane));
		}
		Matrix<Integer> Equations_Ordered(Inequ_on_Ker.nr_of_rows(),dim);
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
		Full_Cone<Integer> Result(ConeDM);
		Result.dual_mode();
		C.set_result(Result);
	}

	*/

	Full_Cone<Integer> FC(ChangeToFullDim.to_sublattice(Matrix<Integer>(Generators)));

	if (computation_type=="triangulation_hilbert_basis") {
		FC.triangulation_hilbert_basis();
	} else if (computation_type=="hilbert_basis") {
		FC.hilbert_basis();
	} else if (computation_type=="support_hyperplanes") {
		FC.support_hyperplanes();
	} else if (computation_type=="support_hyperplanes_pyramid") {
		FC.support_hyperplanes_pyramid();
	} else if (computation_type=="triangulation") {
		FC.support_hyperplanes_triangulation();
	} else if (computation_type=="triangulation_pyramid") {
		FC.support_hyperplanes_triangulation_pyramid();
	} else if (computation_type=="ht1_elements") {
		FC.ht1_elements();
	} else if (computation_type=="hilbert_polynomial") {
		FC.hilbert_polynomial();
	} else if (computation_type=="hilbert_basis_polynomial") {
		FC.hilbert_basis_polynomial();
	} else if (computation_type=="dual") {
		cerr<<"WAHHHH! not implemented yet!";
		throw 1; //TODO implement
	}
}

} //end namespace
