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

#include <list>
#include <math.h>

#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qmap_operations.h"
#include "libQnormaliz/Qconvert.h"
#include "libQnormaliz/Qcone.h"
#include "libQnormaliz/Qfull_cone.h"

// #include "libnormaliz/cone.h"



namespace libQnormaliz {
using namespace std;

// adds the signs inequalities given by Signs to Inequalities
template<typename Number>
Matrix<Number> sign_inequalities(const vector< vector<Number> >& Signs) {
    if (Signs.size() != 1) {
        throw BadInputException("ERROR: Bad signs matrix, has "
                + toString(Signs.size()) + " rows (should be 1)!");
    }
    size_t dim = Signs[0].size();
    Matrix<Number> Inequ(0,dim);
    vector<Number> ineq(dim,0);
    for (size_t i=0; i<dim; i++) {
        Number sign = Signs[0][i];
        if (sign == 1 || sign == -1) {
            ineq[i] = sign;
            Inequ.append(ineq);
            ineq[i] = 0;
        } else if (sign != 0) {
            throw BadInputException("Bad signs matrix, has entry "
                    + toString(sign) + " (should be -1, 1 or 0)!");
        }
    }
    return Inequ;
}

template<typename Number>
Matrix<Number> strict_sign_inequalities(const vector< vector<Number> >& Signs) {
    if (Signs.size() != 1) {
        throw BadInputException("ERROR: Bad signs matrix, has "
                + toString(Signs.size()) + " rows (should be 1)!");
    }
    size_t dim = Signs[0].size();
    Matrix<Number> Inequ(0,dim);
    vector<Number> ineq(dim,0);
    ineq[dim-1]=-1;
    for (size_t i=0; i<dim-1; i++) {    // last component of strict_signs always 0
        Number sign = Signs[0][i];
        if (sign == 1 || sign == -1) {
            ineq[i] = sign;
            Inequ.append(ineq);
            ineq[i] = 0;
        } else if (sign != 0) {
            throw BadInputException("Bad signs matrix, has entry "
                    + toString(sign) + " (should be -1, 1 or 0)!");
        }
    }
    return Inequ;
}

template<typename Number>
vector<vector<Number> > find_input_matrix(const map< InputType, vector< vector<Number> > >& multi_input_data,
                               const InputType type){

    typename map< InputType , vector< vector<Number> > >::const_iterator it;
    it = multi_input_data.find(type);
    if (it != multi_input_data.end())
        return(it->second);

     vector< vector<Number> > dummy;
     return(dummy);
}

template<typename Number>
void insert_column(vector< vector<Number> >& mat, size_t col, Number entry){

    vector<Number> help(mat[0].size()+1);
    for(size_t i=0;i<mat.size();++i){
        for(size_t j=0;j<col;++j)
            help[j]=mat[i][j];
        help[col]=entry;
        for(size_t j=col;j<mat[i].size();++j)
            help[j+1]=mat[i][j];
        mat[i]=help;
    }
}

template<typename Number>
void insert_zero_column(vector< vector<Number> >& mat, size_t col){
    // Number entry=0;
    insert_column<Number>(mat,col,0);
}

template<typename Number>
void Cone<Number>::homogenize_input(map< InputType, vector< vector<Number> > >& multi_input_data){

    typename map< InputType , vector< vector<Number> > >::iterator it;
    it = multi_input_data.begin();
    for(;it!=multi_input_data.end();++it){
        switch(it->first){
            case Type::dehomogenization:
                throw BadInputException("Type dehomogenization not allowed with inhomogeneous input!");
                break;
            case Type::inhom_inequalities: // nothing to do
            case Type::inhom_equations:
            case Type::inhom_congruences:
            case Type::polyhedron:
            case Type::vertices:
            case Type::support_hyperplanes:
            case Type::extreme_rays:
            case Type::open_facets:
            case Type::hilbert_basis_rec_cone:
            case Type::grading:  // already taken care of
                break;
            case Type::strict_inequalities:
                insert_column<Number>(it->second,dim-1,-1);
                break;
            case Type::offset:
                insert_column<Number>(it->second,dim-1,1);
                break;
            default:  // is correct for signs and strict_signs !
                insert_zero_column<Number>(it->second,dim-1);
                break;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
Cone<Number>::Cone(InputType input_type, const vector< vector<Number> >& Input) {
    // convert to a map
    map< InputType, vector< vector<Number> > > multi_input_data;
    multi_input_data[input_type] = Input;
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(InputType type1, const vector< vector<Number> >& Input1,
                    InputType type2, const vector< vector<Number> >& Input2) {
    if (type1 == type2) {
        throw BadInputException("Input types must  pairwise different!");
    }
    // convert to a map
    map< InputType, vector< vector<Number> > > multi_input_data;
    multi_input_data[type1] = Input1;
    multi_input_data[type2] = Input2;
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(InputType type1, const vector< vector<Number> >& Input1,
                    InputType type2, const vector< vector<Number> >& Input2,
                    InputType type3, const vector< vector<Number> >& Input3) {
    if (type1 == type2 || type1 == type3 || type2 == type3) {
        throw BadInputException("Input types must be pairwise different!");
    }
    // convert to a map
    map< InputType, vector< vector<Number> > > multi_input_data;
    multi_input_data[type1] = Input1;
    multi_input_data[type2] = Input2;
    multi_input_data[type3] = Input3;
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(const map< InputType, vector< vector<Number> > >& multi_input_data) {
    process_multi_input(multi_input_data);
}

//---------------------------------------------------------------------------

template<typename Number>
Cone<Number>::Cone(InputType input_type, const Matrix<Number>& Input) {
    // convert to a map
    map< InputType, vector< vector<Number> > >multi_input_data;
    multi_input_data[input_type] = Input.get_elements();
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(InputType type1, const Matrix<Number>& Input1,
                    InputType type2, const Matrix<Number>& Input2) {
    if (type1 == type2) {
        throw BadInputException("Input types must  pairwise different!");
    }
    // convert to a map
    map< InputType, vector< vector<Number> > > multi_input_data;
    multi_input_data[type1] = Input1.get_elements();
    multi_input_data[type2] = Input2.get_elements();
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(InputType type1, const Matrix<Number>& Input1,
                    InputType type2, const Matrix<Number>& Input2,
                    InputType type3, const Matrix<Number>& Input3) {
    if (type1 == type2 || type1 == type3 || type2 == type3) {
        throw BadInputException("Input types must be pairwise different!");
    }
    // convert to a map
    map< InputType, vector< vector<Number> > > multi_input_data;
    multi_input_data[type1] = Input1.get_elements();
    multi_input_data[type2] = Input2.get_elements();
    multi_input_data[type3] = Input3.get_elements();
    process_multi_input(multi_input_data);
}

template<typename Number>
Cone<Number>::Cone(const map< InputType, Matrix<Number> >& multi_input_data_Matrix){
    map< InputType, vector< vector<Number> > > multi_input_data;
    auto it = multi_input_data_Matrix.begin();
    for(; it != multi_input_data_Matrix.end(); ++it){
        multi_input_data[it->first]=it->second.get_elements();
    }
    process_multi_input(multi_input_data);
}

//---------------------------------------------------------------------------


template<typename Number>
void Cone<Number>::process_multi_input(const map< InputType, vector< vector<Number> > >& multi_input_data_const) {
    initialize();
    map< InputType, vector< vector<Number> > > multi_input_data(multi_input_data_const);
    // find basic input type
    bool lattice_ideal_input=false;
    bool inhom_input=false;
    size_t nr_latt_gen=0, nr_cone_gen=0;
    
    auto it = multi_input_data.begin();
    /* for(; it != multi_input_data.end(); ++it)
        for(size_t i=0;i < it->second.size();++i){
            for(size_t j=0;j<it->second[i].size();++j)
                it->second[i][j].canonicalize();
            v_simplify(it->second[i]);
    }*/
    
    inequalities_present=false; //control choice of positive orthant
    
    if(    exists_element(multi_input_data,Type::lattice)
        || exists_element(multi_input_data,Type::lattice_ideal)
        || exists_element(multi_input_data,Type::cone_and_lattice)
        || exists_element(multi_input_data,Type::congruences)
        || exists_element(multi_input_data,Type::inhom_congruences)
        // || exists_element(multi_input_data,Type::dehomogenization)
        || exists_element(multi_input_data,Type::offset)
        || exists_element(multi_input_data,Type::excluded_faces)
        || exists_element(multi_input_data,Type::open_facets)
        || exists_element(multi_input_data,Type::hilbert_basis_rec_cone)
        || exists_element(multi_input_data,Type::strict_inequalities)
        )
        throw BadInputException("Input type not allowed for field coefficients");    

    // NEW: Empty matrix have syntactical influence
    it = multi_input_data.begin();
    for(; it != multi_input_data.end(); ++it) {
        switch (it->first) {
            case Type::inhom_inequalities:
            case Type::inhom_equations:
            case Type::inhom_congruences:
            case Type::strict_inequalities:
            case Type::strict_signs:
                inhom_input=true;
            case Type::signs:
            case Type::inequalities:
            case Type::equations:
            case Type::congruences:
                break;
            case Type::lattice_ideal:
                lattice_ideal_input=true;
                break;
            case Type::polyhedron:
                inhom_input=true;
            case Type::integral_closure:
            case Type::rees_algebra:
            case Type::polytope:
            case Type::cone:
            case Type::subspace:
                nr_cone_gen++;
                break;
            case Type::normalization:
            case Type::cone_and_lattice:
                nr_cone_gen++;
            case Type::lattice:
            case Type::saturation:
                nr_latt_gen++;
                break;
            case Type::vertices:
            case Type::offset:
                inhom_input=true;
            default:
                break;
        }

        switch (it->first) {  // chceck existence of inrqualities
            case Type::inhom_inequalities:
            case Type::strict_inequalities:
            case Type::strict_signs:
            case Type::signs:
            case Type::inequalities:
            case Type::excluded_faces:
                inequalities_present=true;
            default:
                break;
        }

    }

    bool gen_error=false;
    if(nr_cone_gen>2)
        gen_error=true;

    if(nr_cone_gen==2 && (!exists_element(multi_input_data,Type::subspace)
                      || !(exists_element(multi_input_data,Type::cone)
                          || exists_element(multi_input_data,Type::cone_and_lattice)
                          || exists_element(multi_input_data,Type::integral_closure)
                          || exists_element(multi_input_data,Type::normalization) ) )
    )
        gen_error=true;
    
    if(gen_error){
        throw BadInputException("Illegal combination of cone generator types!");
    }
    
    
    if(nr_latt_gen>1){
        throw BadInputException("Only one matrix of lattice generators allowed!");
    }
    if(lattice_ideal_input){
        if(multi_input_data.size() > 2 || (multi_input_data.size()==2 && !exists_element(multi_input_data,Type::grading))){
            throw BadInputException("Only grading allowed with lattice_ideal!");
        }
    }
    if(inhom_input){
        if(exists_element(multi_input_data,Type::dehomogenization) || exists_element(multi_input_data,Type::support_hyperplanes)){
            throw BadInputException("Types dehomogenization and support_hyperplanes not allowed with inhomogeneous input!");
        }
    }
    if(inhom_input || exists_element(multi_input_data,Type::dehomogenization)){
        if(exists_element(multi_input_data,Type::rees_algebra) || exists_element(multi_input_data,Type::polytope)){
            throw BadInputException("Types polytope and rees_algebra not allowed with inhomogeneous input or hehomogenizaion!");
        }
        if(exists_element(multi_input_data,Type::excluded_faces)){
            throw BadInputException("Type excluded_faces not allowed with inhomogeneous input or dehomogenization!");
        }
    }
    if(exists_element(multi_input_data,Type::grading) && exists_element(multi_input_data,Type::polytope)){
           throw BadInputException("No explicit grading allowed with polytope!");
    }
    
    // remove empty matrices
    it = multi_input_data.begin();
    for(; it != multi_input_data.end();) {
        if (it->second.size() == 0)
            multi_input_data.erase(it++);
        else
            ++it;
    }

    if(multi_input_data.size()==0){
        throw BadInputException("All input matrices empty!");
    }

    //determine dimension
    it = multi_input_data.begin();
    size_t inhom_corr = 0; // correction in the inhom_input case
    if (inhom_input) inhom_corr = 1;
    dim = it->second.front().size() - type_nr_columns_correction(it->first) + inhom_corr;

    // We now process input types that are independent of generators, constraints, lattice_ideal
    // check for excluded faces
    ExcludedFaces = find_input_matrix(multi_input_data,Type::excluded_faces);
    
    // PreComputedSupportHyperplanes = find_input_matrix(multi_input_data,Type::support_hyperplanes);
    
    // check for a grading
    vector< vector<Number> > lf = find_input_matrix(multi_input_data,Type::grading);
    if (lf.size() > 1) {
        throw BadInputException("Bad grading, has "
                + toString(lf.size()) + " rows (should be 1)!");
    }
    if(lf.size()==1){
        if(inhom_input)
            lf[0].push_back(0); // first we extend grading trivially to have the right dimension
        setGrading (lf[0]);     // will eventually be set in full_cone.cpp

    }


    // check consistence of dimension
    it = multi_input_data.begin();
    size_t test_dim;
    for (; it != multi_input_data.end(); ++it) {
        test_dim = it->second.front().size() - type_nr_columns_correction(it->first) + inhom_corr;
        if (test_dim != dim) {
            throw BadInputException("Inconsistent dimensions in input!");
        }
    }

    if(inhom_input)
        homogenize_input(multi_input_data);
    
    // check for dehomogenization
    lf = find_input_matrix(multi_input_data,Type::dehomogenization);
    if (lf.size() > 1) {
        throw BadInputException("Bad dehomogenization, has "
                + toString(lf.size()) + " rows (should be 1)!");
    }
    if(lf.size()==1){
        setDehomogenization(lf[0]);
    }

    // now we can unify implicit and explicit truncation
    // Note: implicit and explicit truncation have already been excluded
    if (inhom_input) {
        Dehomogenization.resize(dim,0),
        Dehomogenization[dim-1]=1;
        is_Computed.set(ConeProperty::Dehomogenization);
    }
    if(isComputed(ConeProperty::Dehomogenization))
        inhomogeneous=true;

    if(lattice_ideal_input){
        prepare_input_lattice_ideal(multi_input_data);
    }

    Matrix<Number> LatticeGenerators(0,dim);
    prepare_input_generators(multi_input_data, LatticeGenerators);

    Matrix<Number> Equations(0,dim), Congruences(0,dim+1);
    Matrix<Number> Inequalities(0,dim);
    prepare_input_constraints(multi_input_data,Equations,Congruences,Inequalities);

    // set default values if necessary
    if(inhom_input && LatticeGenerators.nr_of_rows()!=0 && !exists_element(multi_input_data,Type::offset)){
        vector<Number> offset(dim);
        offset[dim-1]=1;
        LatticeGenerators.append(offset);
    }
    if(inhom_input &&  Generators.nr_of_rows()!=0 && !exists_element(multi_input_data,Type::vertices) 
                && !exists_element(multi_input_data,Type::polyhedron)){
        vector<Number> vertex(dim);
        vertex[dim-1]=1;
        Generators.append(vertex);
    }

    if(Inequalities.nr_of_rows()>0 && Generators.nr_of_rows()>0){ // eliminate superfluous inequalities
        vector<key_t> essential;
        for(size_t i=0;i<Inequalities.nr_of_rows();++i){
            for (size_t j=0;j<Generators.nr_of_rows();++j){
                if(v_scalar_product(Inequalities[i],Generators[j])<0){
                    essential.push_back(i);
                    break;
                }
            }
        }
        if(essential.size()<Inequalities.nr_of_rows())
            Inequalities=Inequalities.submatrix(essential);
    }

    // cout << "Ineq " << Inequalities.nr_of_rows() << endl;

    process_lattice_data(LatticeGenerators,Congruences,Equations);

    bool cone_sat_eq=no_lattice_restriction;
    bool cone_sat_cong=no_lattice_restriction;

    // cout << "nolatrest " << no_lattice_restriction << endl;

    if(Inequalities.nr_of_rows()==0 && Generators.nr_of_rows()!=0){
        if(!no_lattice_restriction){
            cone_sat_eq=true;
            for(size_t i=0;i<Generators.nr_of_rows() && cone_sat_eq;++i)
                for(size_t j=0;j<Equations.nr_of_rows()  && cone_sat_eq ;++j)
                    if(v_scalar_product(Generators[i],Equations[j])!=0){
                        cone_sat_eq=false;
            }
        }

        if(cone_sat_eq && !cone_sat_cong){ // multiply generators by anniullator mod sublattice
            for(size_t i=0;i<Generators.nr_of_rows();++i)
                v_scalar_multiplication(Generators[i],BasisChange.getAnnihilator());
            cone_sat_cong=true;
        }
    }

    if((Inequalities.nr_of_rows()!=0 || !cone_sat_eq) && Generators.nr_of_rows()!=0){
        Sublattice_Representation<Number> ConeLatt(Generators,true);
        Full_Cone<Number> TmpCone(ConeLatt.to_sublattice(Generators));
        TmpCone.dualize_cone();
        Inequalities.append(ConeLatt.from_sublattice_dual(TmpCone.Support_Hyperplanes));
        Generators=Matrix<Number>(0,dim); // Generators now converted into inequalities
    }

    assert(Inequalities.nr_of_rows()==0 || Generators.nr_of_rows()==0);    

    if(Generators.nr_of_rows()==0)
        prepare_input_type_4(Inequalities); // inserts default inequalties if necessary
    else{
        is_Computed.set(ConeProperty::Generators);
        is_Computed.set(ConeProperty::Sublattice); 
    }
    
    checkDehomogenization();
    checkGrading();
    
    setWeights();  // make matrix of weights for sorting
    
    if(exists_element(multi_input_data,Type::support_hyperplanes)){
        // SupportHyperplanes=PreComputedSupportHyperplanes;
        SupportHyperplanes = find_input_matrix(multi_input_data,Type::support_hyperplanes);
        is_Computed.set(ConeProperty::SupportHyperplanes);
    }
    
    if(exists_element(multi_input_data,Type::extreme_rays)){
        // SupportHyperplanes=PreComputedSupportHyperplanes;
        Generators = find_input_matrix(multi_input_data,Type::extreme_rays);
        is_Computed.set(ConeProperty::Generators);
        set_extreme_rays(vector<bool>(Generators.nr_of_rows(),true));
    }
    
    BasisChangePointed=BasisChange;
    
    is_Computed.set(ConeProperty::IsInhomogeneous);
    is_Computed.set(ConeProperty::EmbeddingDim);
    
    if(isComputed(ConeProperty::Generators)){
        vector<Number> Grad;
        if(inhomogeneous)
            Grad=Dehomogenization;
        else
            Grad=Grading;
        for(size_t i=0;i<Generators.nr_of_rows();++i)
            v_simplify(Generators[i],Grad);
    }

    /* if(ExcludedFaces.nr_of_rows()>0){ // Nothing to check anymore
        check_excluded_faces();
    } */

    /*
    cout <<"-----------------------" << endl;
    cout << "Gen " << endl;
    Generators.pretty_print(cout);
    cout << "Supp " << endl;
    SupportHyperplanes.pretty_print(cout);
    cout << "A" << endl;
    BasisChange.get_A().pretty_print(cout);
    cout << "B" << endl;
    BasisChange.get_B().pretty_print(cout);
    cout <<"-----------------------" << endl;
    */
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::setGrading (const vector<Number>& lf) {
    
    if (isComputed(ConeProperty::Grading) && Grading == lf) {
        return;
    }
    
    if (lf.size() != dim) {
        throw BadInputException("Grading linear form has wrong dimension "
                + toString(lf.size()) + " (should be " + toString(dim) + ")");
    }
    
    Grading = lf;
    checkGrading();
}

template<typename Number>
void Cone<Number>::checkGrading () {
    
    if (isComputed(ConeProperty::Grading) || Grading.size()==0) {
        return;
    }
    
    bool positively_graded=true;
    bool nonnegative=true;
    size_t neg_index=0;
    Number neg_value;
    if (Generators.nr_of_rows() > 0) {
        vector<Number> degrees = Generators.MxV(Grading);
        for (size_t i=0; i<degrees.size(); ++i) {
            if (degrees[i]<=0 && (!inhomogeneous || v_scalar_product(Generators[i],Dehomogenization)==0)) { 
                // in the inhomogeneous case: test only generators of tail cone
                positively_graded=false;;
                if(degrees[i]<0){
                    nonnegative=false;
                    neg_index=i;
                    neg_value=degrees[i];
                }
            }
        }
        if(positively_graded){
            //vector<Number> test_grading=BasisChange.to_sublattice_dual_no_div(Grading);
           //  GradingDenom=v_make_prime(test_grading);
            GradingDenom=1;
        }
        else
            GradingDenom = 1; 
    } else {
        GradingDenom = 1;
    }

    if (isComputed(ConeProperty::Generators)){        
        if(!nonnegative){
            throw BadInputException("Grading gives negative value "
                    + toString(neg_value) + " for generator "
                    + toString(neg_index+1) + "!");
        }
        if(positively_graded){
            is_Computed.set(ConeProperty::Grading);
            is_Computed.set(ConeProperty::GradingDenom);            
        }
    }
    
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::prepare_input_constraints(const map< InputType, vector< vector<Number> > >& multi_input_data,
    Matrix<Number>& Equations, Matrix<Number>& Congruences, Matrix<Number>& Inequalities) {

    Matrix<Number> Signs(0,dim), StrictSigns(0,dim);

    SupportHyperplanes=Matrix<Number>(0,dim);

    typename map< InputType , vector< vector<Number> > >::const_iterator it=multi_input_data.begin();

    it = multi_input_data.begin();
    for (; it != multi_input_data.end(); ++it) {

        switch (it->first) {
            case Type::strict_inequalities:
            case Type::inequalities:
            case Type::inhom_inequalities:
            case Type::excluded_faces:
                Inequalities.append(it->second);
                break;
            case Type::equations:
            case Type::inhom_equations:
                Equations.append(it->second);
                break;
            case Type::congruences:
            case Type::inhom_congruences:
                Congruences.append(it->second);
                break;
            case Type::signs:
                Signs = sign_inequalities(it->second);
                break;
            case Type::strict_signs:
                StrictSigns = strict_sign_inequalities(it->second);
                break;
            default:
                break;
        }
    }
    if(!BC_set) compose_basis_change(Sublattice_Representation<Number>(dim));
    Matrix<Number> Help(Signs);  // signs first !!
    Help.append(StrictSigns);   // then strict signs
    Help.append(Inequalities);
    Inequalities=Help;
}

//---------------------------------------------------------------------------
template<typename Number>
void Cone<Number>::prepare_input_generators(map< InputType, vector< vector<Number> > >& multi_input_data, Matrix<Number>& LatticeGenerators) {

    if(exists_element(multi_input_data,Type::vertices)){
        for(size_t i=0;i<multi_input_data[Type::vertices].size();++i)
            if(multi_input_data[Type::vertices][i][dim-1] <= 0) {
                throw BadInputException("Vertex has non-positive denominator!");
            }
    }

    if(exists_element(multi_input_data,Type::polyhedron)){
        for(size_t i=0;i<multi_input_data[Type::polyhedron].size();++i)
            if(multi_input_data[Type::polyhedron][i][dim-1] < 0) {
                throw BadInputException("Polyhedron vertex has negative denominator!");
            }
    }

    typename map< InputType , vector< vector<Number> > >::const_iterator it=multi_input_data.begin();
    // find specific generator type -- there is only one, as checked already

    normalization=false;
    
    // check for subspace
    BasisMaxSubspace = find_input_matrix(multi_input_data,Type::subspace);
    if(BasisMaxSubspace.nr_of_rows()==0)
        BasisMaxSubspace=Matrix<Number>(0,dim);
    
    vector<Number> neg_sum_subspace(dim,0);
    for(size_t i=0;i<BasisMaxSubspace.nr_of_rows();++i)
        neg_sum_subspace=v_add(neg_sum_subspace,BasisMaxSubspace[i]);
    v_scalar_multiplication<Number>(neg_sum_subspace,-1);
    

    Generators=Matrix<Number>(0,dim);
    for(; it != multi_input_data.end(); ++it) {
        switch (it->first) {
            case Type::normalization:
            case Type::cone_and_lattice:
                normalization=true;
                LatticeGenerators.append(it->second);
                if(BasisMaxSubspace.nr_of_rows()>0)
                    LatticeGenerators.append(BasisMaxSubspace);
            case Type::vertices:
            case Type::polyhedron:
            case Type::cone:
            case Type::integral_closure:
                Generators.append(it->second);
                break;
            case Type::subspace:
                Generators.append(it->second);
                Generators.append(neg_sum_subspace);
                break;
            case Type::polytope:
                Generators.append(prepare_input_type_2(it->second));
                break;
            case Type::rees_algebra:
                Generators.append(prepare_input_type_3(it->second));
                break;
            case Type::lattice:
                LatticeGenerators.append(it->second);
                break;
            case Type::saturation:
                LatticeGenerators.append(it->second);
                LatticeGenerators.saturate();
                break;
            case Type::offset:
                if(it->second.size()>1){
                  throw BadInputException("Only one offset allowed!");
                }
                LatticeGenerators.append(it->second);
                break;
            default: break;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::process_lattice_data(const Matrix<Number>& LatticeGenerators, Matrix<Number>& Congruences, Matrix<Number>& Equations) {

    if(!BC_set)
        compose_basis_change(Sublattice_Representation<Number>(dim));

    bool no_constraints=(Congruences.nr_of_rows()==0) && (Equations.nr_of_rows()==0);
    bool only_cone_gen=(Generators.nr_of_rows()!=0) && no_constraints && (LatticeGenerators.nr_of_rows()==0);

    no_lattice_restriction=true;

    if(only_cone_gen){
        Sublattice_Representation<Number> Basis_Change(Generators,true);
        compose_basis_change(Basis_Change);
        return;
    }

    if(normalization && no_constraints){
        Sublattice_Representation<Number> Basis_Change(Generators,false);
        compose_basis_change(Basis_Change);
        return;
    }

    no_lattice_restriction=false;

    if(Generators.nr_of_rows()!=0){
        Equations.append(Generators.kernel());
    }

    if(LatticeGenerators.nr_of_rows()!=0){
        Sublattice_Representation<Number> GenSublattice(LatticeGenerators,false);
        if((Equations.nr_of_rows()==0) && (Congruences.nr_of_rows()==0)){
            compose_basis_change(GenSublattice);
            return;
        }
        // Congruences.append(GenSublattice.getCongruencesMatrix());
        Equations.append(GenSublattice.getEquationsMatrix());
    }

    /* if (Congruences.nr_of_rows() > 0) {
        bool zero_modulus;
        Matrix<Number> Ker_Basis=Congruences.solve_congruences(zero_modulus);
        if(zero_modulus) {
            throw BadInputException("Modulus 0 in congruence!");
        }
        Sublattice_Representation<Number> Basis_Change(Ker_Basis,false);
        compose_basis_change(Basis_Change);
    }*/

    if (Equations.nr_of_rows()>0) {
        Matrix<Number> Ker_Basis=BasisChange.to_sublattice_dual(Equations).kernel();
        Sublattice_Representation<Number> Basis_Change(Ker_Basis,true);
        compose_basis_change(Basis_Change);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::prepare_input_type_4(Matrix<Number>& Inequalities) {

    if (!inequalities_present) {
        if (verbose) {
            verboseOutput() << "No inequalities specified in constraint mode, using non-negative orthant." << endl;
        }
        if(inhomogeneous){
            vector<Number> test(dim);
            test[dim-1]=1;
            size_t matsize=dim;
            if(test==Dehomogenization) // in this case "last coordinate >= 0" will come in through the dehomogenization
                matsize=dim-1;   // we don't check for any other coincidence
            Inequalities= Matrix<Number>(matsize,dim);
            for(size_t j=0;j<matsize;++j)
                Inequalities[j][j]=1;
        }
        else
            Inequalities = Matrix<Number>(dim);
    }
    if(inhomogeneous)
        SupportHyperplanes.append(Dehomogenization);
    SupportHyperplanes.append(Inequalities);
}


//---------------------------------------------------------------------------

/* polytope input */
template<typename Number>
Matrix<Number> Cone<Number>::prepare_input_type_2(const vector< vector<Number> >& Input) {
    size_t j;
    size_t nr = Input.size();
    //append a column of 1
    Matrix<Number> Generators(nr, dim);
    for (size_t i=0; i<nr; i++) {
        for (j=0; j<dim-1; j++)
            Generators[i][j] = Input[i][j];
        Generators[i][dim-1]=1;
    }
    // use the added last component as grading
    Grading = vector<Number>(dim,0);
    Grading[dim-1] = 1;
    is_Computed.set(ConeProperty::Grading);
    GradingDenom=1;
    is_Computed.set(ConeProperty::GradingDenom);
    return Generators;
}

//---------------------------------------------------------------------------

/* rees input */
template<typename Number>
Matrix<Number> Cone<Number>::prepare_input_type_3(const vector< vector<Number> >& InputV) {
    Matrix<Number> Input(InputV);
    int i,j,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
    // create cone generator matrix
    Matrix<Number> Full_Cone_Generators(nr_rows+nr_columns,nr_columns+1,0);
    for (i = 0; i < nr_columns; i++) {
        Full_Cone_Generators[i][i]=1;
    }
    for(i=0; i<nr_rows; i++){
        Full_Cone_Generators[i+nr_columns][nr_columns]=1;
        for(j=0; j<nr_columns; j++) {
            Full_Cone_Generators[i+nr_columns][j]=Input[i][j];
        }
    }

    return Full_Cone_Generators;
}


//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::prepare_input_lattice_ideal(map< InputType, vector< vector<Number> > >& multi_input_data) {

    Matrix<Number> Binomials(find_input_matrix(multi_input_data,Type::lattice_ideal));

    if (Grading.size()>0) {
        //check if binomials are homogeneous
        vector<Number> degrees = Binomials.MxV(Grading);
        for (size_t i=0; i<degrees.size(); ++i) {
            if (degrees[i]!=0) {
                throw BadInputException("Grading gives non-zero value "
                        + toString(degrees[i]) + " for binomial "
                        + toString(i+1) + "!");
            }
            if (Grading[i] <0) {
                throw BadInputException("Grading gives negative value "
                        + toString(Grading[i]) + " for generator "
                        + toString(i+1) + "!");
            }
        }
    }

    Matrix<Number> Gens=Binomials.kernel().transpose();
    Full_Cone<Number> FC(Gens);
    FC.verbose=verbose;
    if (verbose) verboseOutput() << "Computing a positive embedding..." << endl;

    FC.dualize_cone();
    Matrix<Number> Supp_Hyp=FC.getSupportHyperplanes().sort_lex();
    Matrix<Number> Selected_Supp_Hyp_Trans=(Supp_Hyp.submatrix(Supp_Hyp.max_rank_submatrix_lex())).transpose();
    Matrix<Number> Positive_Embedded_Generators=Gens.multiplication(Selected_Supp_Hyp_Trans);
    // GeneratorsOfToricRing = Positive_Embedded_Generators;
    // is_Computed.set(ConeProperty::GeneratorsOfToricRing);
    dim = Positive_Embedded_Generators.nr_of_columns();
    multi_input_data.insert(make_pair(Type::normalization,Positive_Embedded_Generators.get_elements())); // this is the cone defined by the binomials

    if (Grading.size()>0) {
        // solve GeneratorsOfToricRing * grading = old_grading
        Number dummyDenom;
        // Grading must be set directly since map entry has been processed already
        Grading = Positive_Embedded_Generators.solve_rectangular(Grading,dummyDenom);
        if (Grading.size() != dim) {
            errorOutput() << "Grading could not be transferred!"<<endl;
            is_Computed.set(ConeProperty::Grading, false);
        }
    }
}

/* only used by the constructors */
template<typename Number>
void Cone<Number>::initialize() {
    BC_set=false;
    is_Computed = bitset<ConeProperty::EnumSize>();  //initialized to false
    dim = 0;
    unit_group_index = 1;
    inhomogeneous=false;
    triangulation_is_nested = false;
    triangulation_is_partial = false;
    verbose = libQnormaliz::verbose; //take the global default
    dual_original_generators=false;
    
    set_parallelization();
    
    /* if (using_GMP<Number>()) {
        change_integer_type = true;
    } else { */
        change_integer_type = false;
    // }
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::compose_basis_change(const Sublattice_Representation<Number>& BC) {
    if (BC_set) {
        BasisChange.compose(BC);
    } else {
        BasisChange = BC;
        BC_set = true;
    }
}
//---------------------------------------------------------------------------
template<typename Number>
void Cone<Number>::check_precomputed_support_hyperplanes(){

    if (isComputed(ConeProperty::Generators)) {
        // check if the inequalities are at least valid
        // if (PreComputedSupportHyperplanes.nr_of_rows() != 0) {
            Number sp;
            for (size_t i = 0; i < Generators.nr_of_rows(); ++i) {
                for (size_t j = 0; j < PreComputedSupportHyperplanes.nr_of_rows(); ++j) {
                    if ((sp = v_scalar_product(Generators[i], PreComputedSupportHyperplanes[j])) < 0) {
                        throw BadInputException("Precomputed inequality " + toString(j)
                                + " is not valid for generator " + toString(i)
                                + " (value " + toString(sp) + ")");
                    }
                }
            }
        // }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
bool Cone<Number>::setVerbose (bool v) {
    //we want to return the old value
    bool old = verbose;
    verbose = v;
    return old;
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::checkDehomogenization () {
    if(Dehomogenization.size()>0){
        vector<Number> test=Generators.MxV(Dehomogenization);
        for(size_t i=0;i<test.size();++i)
            if(test[i]<0){
                throw BadInputException(
                        "Dehomogenization has has negative value on generator "
                        + toString(Generators[i]));
            }
    }
}


//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::setWeights () {

    if(WeightsGrad.nr_of_columns()!=dim){
        WeightsGrad=Matrix<Number> (0,dim);  // weight matrix for ordering
    }
    if(Grading.size()>0 && WeightsGrad.nr_of_rows()==0)
        WeightsGrad.append(Grading);
    GradAbs=vector<bool>(WeightsGrad.nr_of_rows(),false);
}
//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::setDehomogenization (const vector<Number>& lf) {
    if (lf.size() != dim) {
        throw BadInputException("Dehomogenizing linear form has wrong dimension "
                + toString(lf.size()) + " (should be " + toString(dim) + ")");
    }
    Dehomogenization=lf;
    is_Computed.set(ConeProperty::Dehomogenization);
}

//---------------------------------------------------------------------------

/* check what is computed */

/*template<typename Number>
const renf_class* Cone<Number>::getRenf() const {
    return Renf;
}*/


template<typename Number>
bool Cone<Number>::isComputed(ConeProperty::Enum prop) const {
    return is_Computed.test(prop);
}

template<typename Number>
bool Cone<Number>::isComputed(ConeProperties CheckComputed) const {
    return CheckComputed.reset(is_Computed).any();
}


/* getter */

template<typename Number>
size_t Cone<Number>::getRank() {
    compute(ConeProperty::Sublattice);
    return BasisChange.getRank();
}


template<typename Number>
size_t Cone<Number>::getRecessionRank() {
    compute(ConeProperty::RecessionRank);
    return recession_rank;
}

template<typename Number>
long Cone<Number>::getAffineDim() {
    compute(ConeProperty::AffineDim);
    return affine_dim;
}

template<typename Number>
const Sublattice_Representation<Number>& Cone<Number>::getSublattice() {
    compute(ConeProperty::Sublattice);
    return BasisChange;
}

template<typename Number>
const vector< vector<Number> >& Cone<Number>::getMaximalSubspace() {
    compute(ConeProperty::MaximalSubspace);
    return BasisMaxSubspace.get_elements();
}
template<typename Number>
const Matrix<Number>& Cone<Number>::getMaximalSubspaceMatrix() {
    compute(ConeProperty::MaximalSubspace);
    return BasisMaxSubspace;
}
template<typename Number>
size_t Cone<Number>::getDimMaximalSubspace() {
    compute(ConeProperty::MaximalSubspace);
    return BasisMaxSubspace.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getGeneratorsMatrix() {
    compute(ConeProperty::Generators);
    return Generators;
}

template<typename Number>
const vector< vector<Number> >& Cone<Number>::getGenerators() {
    compute(ConeProperty::Generators);
    return Generators.get_elements();
}

template<typename Number>
size_t Cone<Number>::getNrGenerators() {
    compute(ConeProperty::Generators);
    return Generators.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getExtremeRaysMatrix() {
    compute(ConeProperty::ExtremeRays);
    return ExtremeRays;
}
template<typename Number>
const vector< vector<Number> >& Cone<Number>::getExtremeRays() {
    compute(ConeProperty::ExtremeRays);
    return ExtremeRays.get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrExtremeRays() {
    compute(ConeProperty::ExtremeRays);
    return ExtremeRays.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getVerticesOfPolyhedronMatrix() {
    compute(ConeProperty::VerticesOfPolyhedron);
    return VerticesOfPolyhedron;
}
template<typename Number>
const vector< vector<Number> >& Cone<Number>::getVerticesOfPolyhedron() {
    compute(ConeProperty::VerticesOfPolyhedron);
    return VerticesOfPolyhedron.get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrVerticesOfPolyhedron() {
    compute(ConeProperty::VerticesOfPolyhedron);
    return VerticesOfPolyhedron.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getSupportHyperplanesMatrix() {
    compute(ConeProperty::SupportHyperplanes);
    return SupportHyperplanes;
}
template<typename Number>
const vector< vector<Number> >& Cone<Number>::getSupportHyperplanes() {
    compute(ConeProperty::SupportHyperplanes);
    return SupportHyperplanes.get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrSupportHyperplanes() {
    compute(ConeProperty::SupportHyperplanes);
    return SupportHyperplanes.nr_of_rows();
}

template<typename Number>
map< InputType , vector< vector<Number> > > Cone<Number>::getConstraints () {
    compute(ConeProperty::Sublattice, ConeProperty::SupportHyperplanes);
    map<InputType, vector< vector<Number> > > c;
    c[Type::inequalities] = SupportHyperplanes.get_elements();
    c[Type::equations] = BasisChange.getEquations();
    // c[Type::congruences] = BasisChange.getCongruences();
    return c;
}

template<typename Number>
const vector< pair<vector<key_t>,Number> >& Cone<Number>::getTriangulation() {
    compute(ConeProperty::Triangulation);
    return Triangulation;
}

template<typename Number>
const vector<vector<bool> >& Cone<Number>::getOpenFacets() {
    compute(ConeProperty::ConeDecomposition);
    return OpenFacets;
}

template<typename Number>
size_t Cone<Number>::getTriangulationSize() {
    compute(ConeProperty::TriangulationSize);
    return TriangulationSize;
}

template<typename Number>
Number Cone<Number>::getTriangulationDetSum() {
    compute(ConeProperty::TriangulationDetSum);
    return TriangulationDetSum;
}

template<typename Number>
vector<Number> Cone<Number>::getDehomogenization() {
    compute(ConeProperty::Dehomogenization);
    return Dehomogenization;
}

template<typename Number>
bool Cone<Number>::isPointed() {
    compute(ConeProperty::IsPointed);
    return pointed;
}

template<typename Number>
bool Cone<Number>::isInhomogeneous() {
    return inhomogeneous;
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getModuleGeneratorsMatrix() {
    compute(ConeProperty::ModuleGenerators);
    return ModuleGenerators;
}
template<typename Number>
const vector< vector<Number> >& Cone<Number>::getModuleGenerators() {
    compute(ConeProperty::ModuleGenerators);
    return ModuleGenerators.get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrModuleGenerators() {
    compute(ConeProperty::ModuleGenerators);
    return ModuleGenerators.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getDeg1ElementsMatrix() {
    compute(ConeProperty::Deg1Elements);
    return Deg1Elements;
}
template<typename Number>
const vector< vector<Number> >& Cone<Number>::getDeg1Elements() {
    compute(ConeProperty::Deg1Elements);
    return Deg1Elements.get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrDeg1Elements() {
    compute(ConeProperty::Deg1Elements);
    return Deg1Elements.nr_of_rows();
}

template<typename Number>
const Matrix<Number>& Cone<Number>::getLatticePointsMatrix() {
    compute(ConeProperty::LatticePoints);
    if(!inhomogeneous)
        return Deg1Elements;
    else
        return ModuleGenerators;        
}

template<typename Number>
const vector< vector<Number> >& Cone<Number>::getLatticePoints() {
    compute(ConeProperty::LatticePoints);
    return getLatticePointsMatrix().get_elements();
}
template<typename Number>
size_t Cone<Number>::getNrLatticePoints() {
    compute(ConeProperty::LatticePoints);
    return getLatticePointsMatrix().nr_of_rows();
}

template<typename Number>
Number Cone<Number>::getVolume() {
    compute(ConeProperty::Volume);
    return volume;
}

template<typename Number>
double Cone<Number>::getEuclideanVolume() {
    compute(ConeProperty::Volume);
    return euclidean_volume;
}

template<typename Number>
Cone<Number>& Cone<Number>::getIntegerHullCone() const {
    return *IntHullCone;
}


// the information about the triangulation will just be returned
// if no triangulation was computed so far they return false
template<typename Number>
bool Cone<Number>::isTriangulationNested() {
    return triangulation_is_nested;
}
template<typename Number>
bool Cone<Number>::isTriangulationPartial() {
    return triangulation_is_partial;
}

//---------------------------------------------------------------------------

template<typename Number>
ConeProperties Cone<Number>::compute(ConeProperty::Enum cp) {
    if (isComputed(cp)) return ConeProperties();
    return compute(ConeProperties(cp));
}

template<typename Number>
ConeProperties Cone<Number>::compute(ConeProperty::Enum cp1, ConeProperty::Enum cp2) {
    return compute(ConeProperties(cp1,cp2));
}

template<typename Number>
ConeProperties Cone<Number>::compute(ConeProperty::Enum cp1, ConeProperty::Enum cp2,
                                      ConeProperty::Enum cp3) {
    return compute(ConeProperties(cp1,cp2,cp3));
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::set_implicit_dual_mode(ConeProperties& ToCompute) {
    
    if(ToCompute.test(ConeProperty::DualMode) || ToCompute.test(ConeProperty::PrimalMode)
                    || ToCompute.test(ConeProperty::ModuleGeneratorsOverOriginalMonoid)
                    || Generators.nr_of_rows()>0 || SupportHyperplanes.nr_of_rows() > 2*dim
                    || SupportHyperplanes.nr_of_rows() 
                            <= BasisChangePointed.getRank()+ 50/(BasisChangePointed.getRank()+1))
        return;
    if(ToCompute.test(ConeProperty::HilbertBasis))
        ToCompute.set(ConeProperty::DualMode);
    if(ToCompute.test(ConeProperty::Deg1Elements) 
            && !(ToCompute.test(ConeProperty::HilbertSeries) || ToCompute.test(ConeProperty::Multiplicity)))
        ToCompute.set(ConeProperty::DualMode);
    return;
}

//---------------------------------------------------------------------------

template<typename Number>
ConeProperties Cone<Number>::compute(ConeProperties ToCompute) {
    
    set_parallelization();
    
    if(ToCompute.test(ConeProperty::GradingIsPositive)){
        if(Grading.size()==0)
            throw BadInputException("No grading declared that could be positive.");
        else
            is_Computed.set(ConeProperty::Grading);       
    }
    
    if(ToCompute.test(ConeProperty::DefaultMode))
        ToCompute.set(ConeProperty::SupportHyperplanes);
    
    change_integer_type=false;
    
    if(BasisMaxSubspace.nr_of_rows()>0 && !isComputed(ConeProperty::MaximalSubspace)){
        BasisMaxSubspace=Matrix<Number>(0,dim);
        compute(ConeProperty::MaximalSubspace);      
    }
    
    
    ToCompute.reset(is_Computed);
    ToCompute.set_preconditions(inhomogeneous);
    
    ToCompute.check_Q_permissible(); // after implications!
    
    ToCompute.prepare_compute_options(inhomogeneous);
    ToCompute.check_sanity(inhomogeneous);

    /* preparation: get generators if necessary */
    compute_generators(ToCompute);

    if (!isComputed(ConeProperty::Generators)) {
        throw FatalException("Could not get Generators.");
    }

    ToCompute.reset(is_Computed); // already computed
    if (ToCompute.none()) {
        return ToCompute;
    }
    
    prepare_volume_computation(ToCompute);

    // the actual computation
    
    if(isComputed(ConeProperty::SupportHyperplanes))
        ToCompute.reset(ConeProperty::DefaultMode);

    if (ToCompute.any()) {
        compute_inner<Number>(ToCompute);
    }
    
    compute_lattice_points_in_polytope(ToCompute);
    
    compute_integer_hull(ToCompute);
    
    complete_sublattice_comp(ToCompute);

    /* check if everything is computed */
    ToCompute.reset(is_Computed); //remove what is now computed
    
    /* if (ToCompute.test(ConeProperty::Deg1Elements) && isComputed(ConeProperty::Grading)) {
        // this can happen when we were looking for a witness earlier
        compute(ToCompute);
    }*/
    
    if (!ToCompute.test(ConeProperty::DefaultMode) && ToCompute.goals().any()) {
        throw NotComputableException(ToCompute.goals());
    }
    ToCompute.reset_compute_options();
    return ToCompute;
}

template<typename Number>
template<typename NumberFC>
void Cone<Number>::compute_inner(ConeProperties& ToCompute) {
    
    if(ToCompute.test(ConeProperty::IsPointed) && Grading.size()==0){
        if (verbose) {
            verboseOutput()<<  "Checking pointedness first"<< endl;
        }
        ConeProperties Dualize;
        Dualize.set(ConeProperty::SupportHyperplanes);
        Dualize.set(ConeProperty::ExtremeRays);
        compute(Dualize);
    }
    
    Matrix<NumberFC> FC_Gens;

    BasisChangePointed.convert_to_sublattice(FC_Gens, Generators);
    Full_Cone<NumberFC> FC(FC_Gens,!ToCompute.test(ConeProperty::ModuleGeneratorsOverOriginalMonoid));
    // !ToCompute.test(ConeProperty::ModuleGeneratorsOverOriginalMonoid) blocks make_prime in full_cone.cpp

    /* activate bools in FC */

    FC.verbose=verbose;

    FC.inhomogeneous=inhomogeneous;

    if (ToCompute.test(ConeProperty::Triangulation)) {
        FC.keep_triangulation = true;
    }
    
    if (ToCompute.test(ConeProperty::Volume)) {
        FC.do_multiplicity= true;
    }
    
    if (ToCompute.test(ConeProperty::ConeDecomposition)) {
        FC.do_cone_dec = true;
    }

    if (ToCompute.test(ConeProperty::TriangulationDetSum) ) {
        FC.do_determinants = true;
    }
    if (ToCompute.test(ConeProperty::TriangulationSize)) {
        FC.do_triangulation = true;
    }
    if (ToCompute.test(ConeProperty::KeepOrder)) {
        FC.keep_order = true;
    }
    
    /* Give extra data to FC */
    if ( isComputed(ConeProperty::ExtremeRays) ) {
        FC.Extreme_Rays_Ind = ExtremeRaysIndicator;
        FC.is_Computed.set(ConeProperty::ExtremeRays);
    }

    if (inhomogeneous){
        BasisChangePointed.convert_to_sublattice_dual_no_div(FC.Truncation, Dehomogenization);
    }

    if (SupportHyperplanes.nr_of_rows()!=0) {
        BasisChangePointed.convert_to_sublattice_dual(FC.Support_Hyperplanes, SupportHyperplanes);
   }
    if (isComputed(ConeProperty::SupportHyperplanes)){
        FC.is_Computed.set(ConeProperty::SupportHyperplanes);
        FC.do_all_hyperplanes = false;
    }
    
    if(isComputed(ConeProperty::Grading)){
        BasisChangePointed.convert_to_sublattice_dual(FC.Grading,Grading);
            FC.is_Computed.set(ConeProperty::Grading);
    }

    /* do the computation */
    
    try {     
        try {
            FC.compute();
        } catch (const NotIntegrallyClosedException& ) {
        }
        is_Computed.set(ConeProperty::Sublattice);
        // make sure we minimize the excluded faces if requested

        extract_data(FC);
        if(isComputed(ConeProperty::IsPointed) && pointed)
            is_Computed.set(ConeProperty::MaximalSubspace);
    } catch(const NonpointedException& ) {
        is_Computed.set(ConeProperty::Sublattice);
        extract_data(FC);
        if(ToCompute.test(ConeProperty::Deg1Elements) || ToCompute.test(ConeProperty::ModuleGenerators)
            || ToCompute.test(ConeProperty::Volume))
            throw NotComputableException("Qnormaliz requuires ointedness for lattice points or volume");
          
        if(verbose){
            verboseOutput() << "Cone not pointed. Restarting computation." << endl;
        }
        FC=Full_Cone<NumberFC>(Matrix<NumberFC>(1)); // to kill the old FC (almost)
        Matrix<Number> Dual_Gen;
        Dual_Gen=BasisChangePointed.to_sublattice_dual(SupportHyperplanes);
        Sublattice_Representation<Number> Pointed(Dual_Gen,true); // sublattice of the dual lattice
        BasisMaxSubspace = BasisChangePointed.from_sublattice(Pointed.getEquationsMatrix());
        BasisMaxSubspace.simplify_rows();
        // check_vanishing_of_grading_and_dehom();
        BasisChangePointed.compose_dual(Pointed);
        is_Computed.set(ConeProperty::MaximalSubspace);        
        // now we get the basis of the maximal subspace
        pointed = (BasisMaxSubspace.nr_of_rows() == 0);
        is_Computed.set(ConeProperty::IsPointed);
        compute_inner<NumberFC>(ToCompute);           
    }
}


template<typename Number>
void Cone<Number>::compute_generators(ConeProperties& ToCompute) {
    //create Generators from SupportHyperplanes
    if (!isComputed(ConeProperty::Generators) && (SupportHyperplanes.nr_of_rows()!=0 ||inhomogeneous)) {
        if (verbose) {
            verboseOutput() << "Computing extreme rays as support hyperplanes of the dual cone:" << endl;
        }

            compute_generators_inner<Number>(ToCompute);

    }
    assert(isComputed(ConeProperty::Generators));
}

template<typename Number>
template<typename NumberFC>
void Cone<Number>::compute_generators_inner(ConeProperties& ToCompute) {
    
    Matrix<Number> Dual_Gen;
    Dual_Gen=BasisChangePointed.to_sublattice_dual(SupportHyperplanes);
    // first we take the quotient of the efficient sublattice modulo the maximal subspace
    Sublattice_Representation<Number> Pointed(Dual_Gen,true); // sublattice of the dual space

    // now we get the basis of the maximal subspace
    if(!isComputed(ConeProperty::MaximalSubspace)){
        BasisMaxSubspace = BasisChangePointed.from_sublattice(Pointed.getEquationsMatrix());
        BasisMaxSubspace.simplify_rows();
        // check_vanishing_of_grading_and_dehom();
        is_Computed.set(ConeProperty::MaximalSubspace);
    }
    if(!isComputed(ConeProperty::IsPointed)){
        pointed = (BasisMaxSubspace.nr_of_rows() == 0);
        is_Computed.set(ConeProperty::IsPointed);
    }
    BasisChangePointed.compose_dual(Pointed); // primal cone now pointed, may not yet be full dimensional

    // restrict the supphyps to efficient sublattice and push to quotient mod subspace
    Matrix<NumberFC> Dual_Gen_Pointed;
    BasisChangePointed.convert_to_sublattice_dual(Dual_Gen_Pointed, SupportHyperplanes);    
    Full_Cone<NumberFC> Dual_Cone(Dual_Gen_Pointed);
    Dual_Cone.verbose=verbose;
    Dual_Cone.do_extreme_rays=true; // we try to find them, need not exist
    if(ToCompute.test(ConeProperty::KeepOrder))
        Dual_Cone.keep_order=true;
    try {     
        Dual_Cone.dualize_cone();
    } catch(const NonpointedException& ){}; // we don't mind if the dual cone is not pointed
    
    if (Dual_Cone.isComputed(ConeProperty::SupportHyperplanes)) {
        //get the extreme rays of the primal cone
        BasisChangePointed.convert_from_sublattice(Generators,
                          Dual_Cone.getSupportHyperplanes());
        is_Computed.set(ConeProperty::Generators);
        
        //get minmal set of support_hyperplanes if possible
        if (Dual_Cone.isComputed(ConeProperty::ExtremeRays)) {            
            Matrix<NumberFC> Supp_Hyp = Dual_Cone.getGenerators().submatrix(Dual_Cone.getExtremeRays());
            BasisChangePointed.convert_from_sublattice_dual(SupportHyperplanes, Supp_Hyp);
            norm_dehomogenization(BasisChangePointed.getRank());
            SupportHyperplanes.sort_lex();
            is_Computed.set(ConeProperty::SupportHyperplanes);
        }
        
        // now the final transformations
        // only necessary if the basis changes computed so far do not make the cone full-dimensional
        // this is equaivalent to the dual cone bot being pointed
        if(!(Dual_Cone.isComputed(ConeProperty::IsPointed) && Dual_Cone.isPointed())){
            // first to full-dimensional pointed
            Matrix<Number> Help;
            Help=BasisChangePointed.to_sublattice(Generators); // sublattice of the primal space
            Sublattice_Representation<Number> PointedHelp(Help,true);
            BasisChangePointed.compose(PointedHelp);
            // second to efficient sublattice
            if(BasisMaxSubspace.nr_of_rows()==0){  // primal cone is pointed and we can copy
                BasisChange=BasisChangePointed;
            }
            else{
                Help=BasisChange.to_sublattice(Generators);
                Help.append(BasisChange.to_sublattice(BasisMaxSubspace));
                Sublattice_Representation<Number> EmbHelp(Help,true); // sublattice of the primal space
                compose_basis_change(EmbHelp);
            }
        }
        is_Computed.set(ConeProperty::Sublattice); // will not be changed anymore
        
        checkGrading();

        setWeights();
        set_extreme_rays(vector<bool>(Generators.nr_of_rows(),true)); // here since they get sorted
        is_Computed.set(ConeProperty::ExtremeRays);
    }
}

template<typename Number>
vector<Sublattice_Representation<Number> > MakeSubAndQuot(const Matrix<Number>& Gen,
                                        const Matrix<Number>& Ker){
    vector<Sublattice_Representation<Number> > Result;                                        
    Matrix<Number> Help=Gen;
    Help.append(Ker);
    Sublattice_Representation<Number> Sub(Help,true);
    Sublattice_Representation<Number> Quot=Sub;
    if(Ker.nr_of_rows()>0){
        Matrix<Number> HelpQuot=Sub.to_sublattice(Ker).kernel();   // kernel here to be interpreted as subspace of the dual
                                                                    // namely the linear forms vanishing on Ker
        Sublattice_Representation<Number> SubToQuot(HelpQuot,true); // sublattice of the dual
        Quot.compose_dual(SubToQuot);
    }
    Result.push_back(Sub);
    Result.push_back(Quot);
    
    return Result;    
}

//---------------------------------------------------------------------------

template<typename Number>
template<typename NumberFC>
void Cone<Number>::extract_data(Full_Cone<NumberFC>& FC) {
    //this function extracts ALL available data from the Full_Cone
    //even if it was in Cone already <- this may change
    //it is possible to delete the data in Full_Cone after extracting it

    if(verbose) {
        verboseOutput() << "transforming data..."<<flush;
    }
    
    if (FC.isComputed(ConeProperty::Generators)) {
        BasisChangePointed.convert_from_sublattice(Generators,FC.getGenerators());
        is_Computed.set(ConeProperty::Generators);
    }
    
    if (FC.isComputed(ConeProperty::IsPointed) && !isComputed(ConeProperty::IsPointed)) {
        pointed = FC.isPointed();
        if(pointed)
            is_Computed.set(ConeProperty::MaximalSubspace);
        is_Computed.set(ConeProperty::IsPointed);
    }    
    

    if (FC.isComputed(ConeProperty::ExtremeRays)) {
        set_extreme_rays(FC.getExtremeRays());
    }
    if (FC.isComputed(ConeProperty::SupportHyperplanes)) {
        extract_supphyps(FC);
        norm_dehomogenization(FC.dim);
        SupportHyperplanes.sort_lex();
        is_Computed.set(ConeProperty::SupportHyperplanes);
    }
    if (FC.isComputed(ConeProperty::TriangulationSize)) {
        TriangulationSize = FC.totalNrSimplices;
        triangulation_is_nested = FC.triangulation_is_nested;
        triangulation_is_partial= FC.triangulation_is_partial;
        is_Computed.set(ConeProperty::TriangulationSize);
        is_Computed.set(ConeProperty::IsTriangulationPartial);
        is_Computed.set(ConeProperty::IsTriangulationNested);
        is_Computed.reset(ConeProperty::Triangulation);
        Triangulation.clear();
    }
    if (FC.isComputed(ConeProperty::TriangulationDetSum)) {
        convert(TriangulationDetSum, FC.detSum);
        is_Computed.set(ConeProperty::TriangulationDetSum);
    }
    
    if (FC.isComputed(ConeProperty::Triangulation)) {
        size_t tri_size = FC.Triangulation.size();
        Triangulation = vector< pair<vector<key_t>, Number> >(tri_size);
        if(FC.isComputed(ConeProperty::ConeDecomposition))
            OpenFacets.resize(tri_size);
        SHORTSIMPLEX<NumberFC> simp;
        for (size_t i = 0; i<tri_size; ++i) {
            simp = FC.Triangulation.front();
            Triangulation[i].first.swap(simp.key);
            // sort(Triangulation[i].first.begin(), Triangulation[i].first.end());
            if (FC.isComputed(ConeProperty::TriangulationDetSum))
                convert(Triangulation[i].second, simp.vol);
            else
                Triangulation[i].second = 0;
            if(FC.isComputed(ConeProperty::ConeDecomposition))
                OpenFacets[i].swap(simp.Excluded);
            FC.Triangulation.pop_front();
        }
        if(FC.isComputed(ConeProperty::ConeDecomposition))
            is_Computed.set(ConeProperty::ConeDecomposition);
        is_Computed.set(ConeProperty::Triangulation);
    }

    if (FC.isComputed(ConeProperty::RecessionRank) && isComputed(ConeProperty::MaximalSubspace)) {
        recession_rank = FC.level0_dim+BasisMaxSubspace.nr_of_rows();
        is_Computed.set(ConeProperty::RecessionRank);
        if (getRank() == recession_rank) {
            affine_dim = -1;
        } else {
            affine_dim = getRank()-1;
        }
        is_Computed.set(ConeProperty::AffineDim);
    }
    
    if(FC.isComputed(ConeProperty::Multiplicity)){
        volume=FC.multiplicity;
        is_Computed.set(ConeProperty::Volume);
        euclidean_volume=approx_to_double(volume);
        for(int i=1;i<dim;++i)
            euclidean_volume/=i;
        euclidean_volume*=euclidean_height;
        
        is_Computed.set(ConeProperty::EuclideanVolume);
    }
    
    /* if (FC.isComputed(ConeProperty::MaximalSubspace) && 
                                   !isComputed(ConeProperty::MaximalSubspace)) {
        BasisChangePointed.convert_from_sublattice(BasisMaxSubspace, FC.Basis_Max_Subspace);
        check_vanishing_of_grading_and_dehom();
        is_Computed.set(ConeProperty::MaximalSubspace);
    }*/

    if (verbose) {
        verboseOutput() << " done." <<endl;
    }
}

//---------------------------------------------------------------------------
template<typename Number>
template<typename NumberFC>
void Cone<Number>::extract_supphyps(Full_Cone<NumberFC>& FC) {
        BasisChangePointed.convert_from_sublattice_dual(SupportHyperplanes, FC.getSupportHyperplanes());
}

template<typename Number>
void Cone<Number>::extract_supphyps(Full_Cone<Number>& FC) {
    if(BasisChangePointed.IsIdentity())
        swap(SupportHyperplanes,FC.Support_Hyperplanes);
    else
        SupportHyperplanes=BasisChangePointed.from_sublattice_dual(FC.getSupportHyperplanes());
}

template<typename Integer>
void Cone<Integer>::norm_dehomogenization(size_t FC_dim){
    if(inhomogeneous && FC_dim<dim){ // make inequality for the inhomogeneous variable appear as dehomogenization
        vector<Integer> dehom_restricted=BasisChangePointed.to_sublattice_dual(Dehomogenization);
        for(size_t i=0;i<SupportHyperplanes.nr_of_rows();++i){
            if(dehom_restricted==BasisChangePointed.to_sublattice_dual(SupportHyperplanes[i])){
                SupportHyperplanes[i]=Dehomogenization;
                break;
            }
        }
    }
}
//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::set_extreme_rays(const vector<bool>& ext) {
    assert(ext.size() == Generators.nr_of_rows());
    ExtremeRaysIndicator=ext;
    vector<bool> choice=ext;
    if (inhomogeneous) {
        // separate extreme rays to rays of the level 0 cone
        // and the verticies of the polyhedron, which are in level >=1
        size_t nr_gen = Generators.nr_of_rows();
        vector<bool> VOP(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            if (ext[i] && v_scalar_product(Generators[i],Dehomogenization) != 0) {
                VOP[i] = true;
                choice[i]=false;
            }
        }
        VerticesOfPolyhedron=Generators.submatrix(VOP);
        VerticesOfPolyhedron.simplify_rows();
        VerticesOfPolyhedron.sort_by_weights(WeightsGrad,GradAbs);
        is_Computed.set(ConeProperty::VerticesOfPolyhedron);
    }
    ExtremeRays=Generators.submatrix(choice);
    ExtremeRays.simplify_rows();
    if(inhomogeneous && !isComputed(ConeProperty::AffineDim) && isComputed(ConeProperty::MaximalSubspace)){
        size_t level0_dim=ExtremeRays.max_rank_submatrix_lex().size();
        recession_rank = level0_dim+BasisMaxSubspace.nr_of_rows();
        is_Computed.set(ConeProperty::RecessionRank);
        if (getRank() == recession_rank) {
            affine_dim = -1;
        } else {
            affine_dim = getRank()-1;
        }
        is_Computed.set(ConeProperty::AffineDim);
        
    }
    ExtremeRays.sort_by_weights(WeightsGrad,GradAbs);
    is_Computed.set(ConeProperty::ExtremeRays);
}

//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::complete_sublattice_comp(ConeProperties& ToCompute) {
    
    if(!isComputed(ConeProperty::Sublattice))
        return;
    is_Computed.set(ConeProperty::Rank);
    if(ToCompute.test(ConeProperty::Equations)){
        BasisChange.getEquationsMatrix(); // just to force computation, ditto below
        is_Computed.set(ConeProperty::Equations);
    }
    /*
    if(ToCompute.test(ConeProperty::Congruences) || ToCompute.test(ConeProperty::ExternalIndex)){
        // BasisChange.getCongruencesMatrix();
        BasisChange.getExternalIndex();
        // is_Computed.set(ConeProperty::Congruences);
        // is_Computed.set(ConeProperty::ExternalIndex);
    }*/
}

//---------------------------------------------------------------------------

/*
// Version with rational approximation and bridge to libnormaliz
// given up in favor of project-and-lift directly for Number, see below.

template<typename Number>
void Cone<Number>::compute_lattice_points_in_polytope(ConeProperties& ToCompute){
    if(isComputed(ConeProperty::ModuleGenerators))
        return;
    if(!ToCompute.test(ConeProperty::ModuleGenerators) && !ToCompute.test(ConeProperty::Deg1Elements))
        return;
    
    if(!isComputed(ConeProperty::Grading) && !isComputed(ConeProperty::Dehomogenization))
        throw BadInputException("Lattice points not computable without grading in the homogeneous case");
        
    compute(ConeProperty::SupportHyperplanes);
    if(!isComputed(ConeProperty::SupportHyperplanes))
        throw FatalException("Could not compute SupportHyperplanes");

    Matrix<Number> Vert;
    if(!inhomogeneous)
       Vert=ExtremeRays;
    else
        Vert=VerticesOfPolyhedron;
    
    vector<Number> LF=Dehomogenization;
    if(!inhomogeneous)
        LF=Grading;
    
    if(!inhomogeneous){
        for(size_t i=0;i<Vert.nr_of_rows();++i)
            if(v_scalar_product(Vert[i],LF) <= 0)
                throw BadInputException("Lattice points not computable for unbounded poyhedra");        
    }
    
    if(getRank()<dim)
        throw BadInputException("Lattice points only computable for full-dimensional poyhedra"); 
    
    vector<mpq_class> ApproxLF=approx_to_mpq(LF);
    for(size_t i=0;i<LF.size();++i)
        if(LF[i]!=ApproxLF[i])
            throw BadInputException("Lattice points only computable with rational dehomogenization or grading");
        
    Matrix<mpq_class> ApproxHyp(SupportHyperplanes.nr_of_rows(),dim); // we make approximations to the support hyperplanes
    for(size_t i=0; i< SupportHyperplanes.nr_of_rows();++i){
        ApproxHyp[i]=approx_to_mpq(SupportHyperplanes[i]);        
    }
    
    for(size_t i=0;i<ApproxHyp.nr_of_rows();++i){ // we modify the approximationsn som that the approximate
        bool not_yet_good;                        // cone conataisn the original vertices
                                                  // by adding small multiples of the grading/dehomogenization ...
        vector<mpq_class> to_add=ApproxLF;
        mpq_class scaled_by=1;
        scaled_by/=100;
        v_scalar_multiplication(to_add, scaled_by); // ... namely 1/100 of it, possibly several times below

        do{                                       
            not_yet_good=false;
            bool first=true;
            for(size_t j=0;j<Vert.nr_of_rows();++j){ 
                Number test=0;
                for(size_t k=0;k<dim;++k)
                    test+=ApproxHyp[i][k]*Vert[j][k];
                if(test<0){
                    not_yet_good=true;
                    ApproxHyp[i]=v_add(ApproxHyp[i],to_add);
                    break;
                }
            }            
        } while(not_yet_good); 
    }
    
    Matrix<mpq_class> LFMat(0,dim); // for the cone constructor
    LFMat.append(ApproxLF);
    
    libnormaliz::Cone<mpz_class> NmzCone(libnormaliz::Type::inequalities,ApproxHyp.get_elements(),
                                         libnormaliz::Type::grading,LFMat.get_elements());
    NmzCone.setVerbose(true);
    NmzCone.compute(libnormaliz::ConeProperty::Deg1Elements, libnormaliz::ConeProperty::Projection);
    vector<vector<mpz_class> > OurDesiredPointsZZ=NmzCone.getDeg1Elements();
    
    vector<vector<Number> > OurDesiredPointsRR; // transfer mpz_class to Number
    for(size_t i=0;i<OurDesiredPointsZZ.size();++i){
        vector<Number> transfer(OurDesiredPointsZZ[i].size());
        for(size_t j=0;j<transfer.size();++j)
            transfer[j]=OurDesiredPointsZZ[i][j];
        OurDesiredPointsRR.push_back(transfer);        
    }
    
    ModuleGenerators=Matrix<Number>(0,dim);
    Deg1Elements=Matrix<Number>(0,dim);
    
    for(size_t i=0;i<OurDesiredPointsRR.size();++i){ // finally discard points outside the original polytope
        bool is_contained=true;
        for(size_t j=0;j<SupportHyperplanes.nr_of_rows();++j){
            if(v_scalar_product(OurDesiredPointsRR[i],SupportHyperplanes[j])<0){
                is_contained=false;
                break;
            }               
        }
        if(is_contained) {
            if(inhomogeneous)
                ModuleGenerators.append(OurDesiredPointsRR[i]);
            else
                Deg1Elements.append(OurDesiredPointsRR[i]);
        }
    }
    
    if(inhomogeneous)    
        is_Computed.set(ConeProperty::ModuleGenerators);
    else
        is_Computed.set(ConeProperty::Deg1Elements);
}

*/

//---------------------------------------------------------------------------


template<typename Number>
void Cone<Number>::compute_lattice_points_in_polytope(ConeProperties& ToCompute){
    if(isComputed(ConeProperty::ModuleGenerators) || isComputed(ConeProperty::Deg1Elements))
        return;
    if(!ToCompute.test(ConeProperty::ModuleGenerators) && !ToCompute.test(ConeProperty::Deg1Elements))
        return;
    
    if(!isComputed(ConeProperty::Grading) && !isComputed(ConeProperty::Dehomogenization))
        throw BadInputException("Lattice points not computable without grading in the homogeneous case");
        
    compute(ConeProperty::SupportHyperplanes);
    if(!isComputed(ConeProperty::SupportHyperplanes))
        throw FatalException("Could not compute SupportHyperplanes");
    
    if(inhomogeneous && ExtremeRays.nr_of_rows()>0 ){
        throw BadInputException("Lattice points not computable for unbounded poyhedra");        
    }
    
    // The same procedure as in cone.cpp, but no approximation, and grading always extra first coordinate
    
    Number MinusOne=-1;
    
    vector<vector<Number> > SuppsHelp=SupportHyperplanes.get_elements();
    Matrix<Number> Equs=BasisChange.getEquationsMatrix();
    for(size_t i=0;i<Equs.nr_of_rows();++i){ // add equations as inequalities
        SuppsHelp.push_back(Equs[i]);
        SuppsHelp.push_back(Equs[i]);
        v_scalar_multiplication(SuppsHelp.back(),MinusOne);
    } 
    Number Zero=0;
    insert_column(SuppsHelp,0,Zero);
    
    // we insert the degree/level into the 0th column
    vector<Number> ExtraEqu(1,-1);
    for(size_t j=0;j<dim;++j){
        if(inhomogeneous)
            ExtraEqu.push_back(Dehomogenization[j]);
        else
            ExtraEqu.push_back(Grading[j]);
    }
    SuppsHelp.push_back(ExtraEqu);
    v_scalar_multiplication(ExtraEqu,MinusOne);
        SuppsHelp.push_back(ExtraEqu);
        
    Matrix<Number> Supps(SuppsHelp);
    
    Matrix<Number> Gens;
    if(inhomogeneous)
        Gens=VerticesOfPolyhedron;
    else
        Gens=ExtremeRays;
    
    Matrix<Number> GradGen(0,dim+1); 
    for(size_t i=0;i<Gens.nr_of_rows();++i){
        vector<Number> gg(dim+1);
        for(size_t j=0;j<dim;++j)
            gg[j+1]=Gens[i][j];
        if(inhomogeneous)
            gg[0]=v_scalar_product(Gens[i],Dehomogenization);
        else
            gg[0]=v_scalar_product(Gens[i],Grading);
        GradGen.append(gg);            
    }
    
    Deg1Elements.resize(0,dim);
    ModuleGenerators.resize(0,dim);
    
    if(inhomogeneous)    
        project_and_lift(ToCompute, GradGen,Supps,false);
    else
        project_and_lift(ToCompute, GradGen,Supps,false);
    
    // In this version, the lattice points are transferresd into the cone
    // in project_and_lift below.
    
    if(inhomogeneous)    
        is_Computed.set(ConeProperty::ModuleGenerators);
    else
        is_Computed.set(ConeProperty::Deg1Elements);
}

//---------------------------------------------------------------------------
template<typename Number>
void Cone<Number>::project_and_lift(ConeProperties& ToCompute, const Matrix<Number>& Gens, Matrix<Number>& Supps, bool float_projection){
    
    // if(verbose)
    //    verboseOutput() << "Starting projection" << endl;
    
    // vector<boost::dynamic_bitset<> > Pair;
   //  vector<boost::dynamic_bitset<> > ParaInPair;
    
    vector< boost::dynamic_bitset<> > Ind;

    //if(!is_parallelotope){
        Ind=vector< boost::dynamic_bitset<> > (Supps.nr_of_rows(), boost::dynamic_bitset<> (Gens.nr_of_rows()));
        for(size_t i=0;i<Supps.nr_of_rows();++i)
            for(size_t j=0;j<Gens.nr_of_rows();++j)
                if(v_scalar_product(Supps[i],Gens[j])==0)
                    Ind[i][j]=true;
    //}
        
    size_t rank=BasisChangePointed.getRank();
    
    Matrix<Number> Verts;
    if(isComputed(ConeProperty::Generators)){
        vector<key_t> choice=identity_key(Gens.nr_of_rows());   //Gens.max_rank_submatrix_lex();
        if(choice.size()>=dim)
            Verts=Gens.submatrix(choice);        
    }
    
    Matrix<mpz_class> Raw(0,Gens.nr_of_columns());

    vector<Number> Dummy;
    // project_and_lift_inner<Number>(Deg1,Supps,Ind,GradingDenom,rank,verbose,true,Dummy);
    ProjectAndLift<Number,mpz_class> PL;
    // if(!is_parallelotope)
        PL=ProjectAndLift<Number,mpz_class>(Supps,Ind,rank);
    //else
     //    PL=ProjectAndLift<Number,Number>(Supps,Pair,ParaInPair,rank);
    PL.set_grading_denom(1);
    PL.set_verbose(verbose);
    PL.set_no_relax(ToCompute.test(ConeProperty::NoRelax));
    PL.set_LLL(false);
    PL.set_vertices(Verts);
    PL.compute();
    PL.put_eg1Points_into(Raw);
    
    for(size_t i=0;i<Raw.nr_of_rows();++i){
        vector<Number> point(dim);
        for(size_t j=0;j<dim;++j){
            point[j]=Raw[i][j+1];            
        }
        if(inhomogeneous)
            ModuleGenerators.append(point);
        else
            Deg1Elements.append(point);
    }
    
    if(verbose)
        verboseOutput() << "Project-and-lift complete" << endl <<
           "------------------------------------------------------------" << endl;
}



//---------------------------------------------------------------------------

template<typename Number>
void Cone<Number>::prepare_volume_computation(ConeProperties& ToCompute){
    if(!ToCompute.test(ConeProperty::Volume))
        return;
    
    if(!inhomogeneous && !isComputed(ConeProperty::Grading))
        throw NotComputableException("Volume neds a grading in the homogeneous case");
    if(getRank()!=dim)
        throw NotComputableException("Qnormaliz rerquires full dimenson for volume");
    vector<Number> Grad;
    if(inhomogeneous)
        Grad=Dehomogenization;
    else
        Grad=Grading;

    /* for(size_t i=0;i<dim;++i)
        if(!Grad[i].is_integer())
            throw NotComputableException("Entries of grading or dehomogenization must be mpzegers for volume");*/
    
    vector<mpz_class> Grad_mpz; //=approx_to_mpq(Grad);
    for(size_t i=0;i<dim;++i)
        Grad_mpz.push_back(Grad[i].get_num());
    for(size_t i=0;i<dim;++i){
        if(Grad[i]!=Grad_mpz[i])
            throw BadInputException("Entries of grading or dehomogenization must be coprime integers for volume");
    }
    //if(libnormaliz::v_make_prime(Grad_mpz)!=1)
    //    throw NotComputableException("Entries of grading or dehomogenization must be coprime integers for volume");
    
    vector<double> Grad_double(dim);
    for(size_t i=0;i<dim;++i)
        // libnormaliz::convert(Grad_double[i],Grad_mpz[i]);
        Grad_double[i]=Grad_mpz[i].get_d();
    
    double norm=v_scalar_product(Grad_double,Grad_double);    
    euclidean_height=sqrt(norm);
}

template<typename Number>
void Cone<Number>::compute_integer_hull(ConeProperties& ToCompute) {
    if(ToCompute.test(ConeProperty::IntegerHull))
        assert(false);
    return;
}

/* template<typename Number>
void Cone<Number>::compute_integer_hull(ConeProperties& ToCompute) {
    
    if(isComputed(ConeProperty::IntegerHull) || !ToCompute.test(ConeProperty::IntegerHull))
        return;
    
    if(verbose){
        verboseOutput() << "Computing integer hull" << endl;
    }

    Matrix<Number> IntHullGen;
    bool IntHullComputable=true;
    size_t nr_extr=0;
    if(inhomogeneous){
        if(!isComputed(ConeProperty::ModuleGenerators))
            IntHullComputable=false;
        IntHullGen=ModuleGenerators;
    }
    else{
        if(!isComputed(ConeProperty::Deg1Elements))
            IntHullComputable=false;
        IntHullGen=Deg1Elements;
    }
    ConeProperties IntHullCompute;
    IntHullCompute.set(ConeProperty::SupportHyperplanes);
    if(!IntHullComputable){
        throw NotComputableException("Integer hull not computable: no integer points available");
    }
    
    if(IntHullGen.nr_of_rows()==0){ // the list of integer points is allowed to be empty
        IntHullGen.append(vector<Number>(dim,0)); // but we need a non-empty input matrix
    }
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    libnormaliz::Matrix<mpz_class> IntHullGen_libnormaliz(IntHullGen.nr_of_rows(),IntHullGen.nr_of_columns());    
    for(size_t i=0;i<IntHullGen.nr_of_rows();++i){
        for(size_t j=0;j<IntHullGen.nr_of_columns();++j){
            IntHullGen_libnormaliz[i][j]=IntHullGen[i][j].get_num();        
        }        
    }
    
    nr_extr=IntHullGen_libnormaliz.extreme_points_first(); // don't need a norm here since all points have degree or level 1

    for(size_t i=0;i<IntHullGen.nr_of_rows();++i){
        for(size_t j=0;j<IntHullGen.nr_of_columns();++j){
            IntHullGen[i][j]=IntHullGen_libnormaliz[i][j];        
        }        
    }

    if(verbose){
        verboseOutput() << nr_extr << " extreme points found"  << endl;
    }
    
    // IntHullGen.pretty_print(cout);
    IntHullCone=new Cone<Number>(Type::cone,IntHullGen.get_elements());
    if(nr_extr!=0)  // we suppress the ordering in full_cone only if we have found few extreme rays
        IntHullCompute.set(ConeProperty::KeepOrder);

    IntHullCone->inhomogeneous=true; // inhomogeneous;
    if(inhomogeneous)
        IntHullCone->Dehomogenization=Dehomogenization;
    else
        IntHullCone->Dehomogenization=Grading;
    IntHullCone->verbose=verbose;
    try{
        IntHullCone->compute(IntHullCompute);
        if(IntHullCone->isComputed(ConeProperty::SupportHyperplanes))
            is_Computed.set(ConeProperty::IntegerHull);
        if(verbose){
            verboseOutput() << "Integer hull finished" << endl;
        }
    }
    catch (const NotComputableException& ){
            errorOutput() << "Error in computation of integer hull" << endl;
    }
}
*/

/*template<typename Number>
void Cone<Number>::set_renf(renf_class *GivenRenf){    
    Renf=GivenRenf;    
}*/

template<typename Integer>
void Cone<Integer>::set_parallelization() {
    
    omp_set_nested(0);
    
    if(thread_limit<0)
        throw BadInputException("Invalid thread limit");
    
    if(parallelization_set){
        if(thread_limit!=0)
            omp_set_num_threads(thread_limit);        
    }
    else{
        if(std::getenv("OMP_NUM_THREADS") == NULL){
            long old=omp_get_max_threads();
            if(old>default_thread_limit)
                set_thread_limit(default_thread_limit);        
            omp_set_num_threads(thread_limit);
        }       
    }
}

template<typename Number>
Cone<Number>::~Cone() {
}



} // end namespace libQnormaliz
