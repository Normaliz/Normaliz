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

//---------------------------------------------------------------------------

#include <stdlib.h>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "Qoutput.h"
#include "libQnormaliz/Qgeneral.h"
#include "libQnormaliz/Qmatrix.h"
#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qmap_operations.h"

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
Output<Number, NumberField>::Output(){
    out=true;
    inv=false;
    ext=false;
    esp=false;
    typ=false;
    egn=false;
    gen=false;
    cst=false;
    tri=false;
    tgn=false;
    ht1=false;
    dec=false;
    lat=false;
    mod=false;
    msp=false;
    lattice_ideal_input = false;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_lattice_ideal_input(bool value){
    lattice_ideal_input=value;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::read() const{
    cout<<"\nname="<<name<<"\n";
    cout<<"\nout="<<out<<"\n";
    cout<<"\ninv="<<inv<<"\n";
    cout<<"\next="<<ext<<"\n";
    cout<<"\nesp="<<esp<<"\n";
    cout<<"\ntyp="<<typ<<"\n";
    cout<<"\negn="<<egn<<"\n";
    cout<<"\ngen="<<gen<<"\n";
    cout<<"\ncst="<<cst<<"\n";
    cout<<"\ntri="<<tri<<"\n";
    cout<<"\ntgn="<<tgn<<"\n";
    cout<<"\nht1="<<ht1<<"\n";
    cout<<"\nResult is:\n";
    Result->print();
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_name(const string& n){
    name=n;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::setCone(Cone<Number> & C) {
    this->Result = &C;
    dim = Result->getEmbeddingDim();
    homogeneous = !Result->isInhomogeneous();
    if (homogeneous) {
        of_cone       = "";
        of_monoid     = "";
        of_polyhedron = "";
    } else {
        of_cone       = " of recession cone";
        of_monoid     = " of recession monoid";
        of_polyhedron = " of polyhedron (homogenized)";
    }
}

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_renf(ostream & os) const{
    
}

#ifdef ENFNORMALIZ
template<>
void Output<renf_elem_class, renf_class>::write_renf(ostream & os) const{
    os << "Real embedded number field:" << endl;
    os << *Renf << endl;  
    
}
#endif

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_renf(NumberField *renf){
    
}

#ifdef ENFNORMALIZ
template<>
void Output<renf_elem_class, renf_class>::set_renf(renf_class *renf){
    
    Renf=renf;
    
}
#endif

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_out(const bool& flag){
    out=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_inv(const bool& flag){
    inv=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_ext(const bool& flag){
    ext=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_esp(const bool& flag){
    esp=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_typ(const bool& flag){
    typ=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_egn(const bool& flag){
    egn=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_gen(const bool& flag){
    gen=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_cst(const bool& flag){
    cst=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_tri(const bool& flag) {
    tri=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_tgn(const bool& flag) {
    tgn=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_ht1(const bool& flag) {
    ht1=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_dec(const bool& flag) {
    dec=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_mod(const bool& flag) {
    mod=flag;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_lat(const bool& flag) {
    lat=flag;
}


//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_msp(const bool& flag) {
    msp=flag;
}
//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_extra_files(){
    out=true;
    inv=true;
    gen=true;
    cst=true;
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::set_write_all_files(){
    out=true;
    inv=true;
    ext=true;
    esp=true;
    typ=true;
    egn=true;
    gen=true;
    cst=true;
    ht1=true;
    lat=true;
    mod=true;
    msp=true;
}



//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_ext(const Matrix<Number>& M) const{
    if (ext==true) {
        M.print(name,"ext");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_mod(const Matrix<Number>& M) const{
    if (mod==true) {
        M.print(name,"mod");
    }
}


//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_lat(const Matrix<Number>& M) const{
    if (ext==true) {
        M.print(name,"lat");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_esp(const Matrix<Number>& M) const{
    if (esp==true) {
        M.print(name,"esp");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_typ(const Matrix<Number>& M) const{
    if (typ==true) {
        M.print(name,"typ");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_egn(const Matrix<Number>& M) const {
    if (egn==true) {
        M.print(name,"egn");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_gen(const Matrix<Number>& M) const {
    if (gen==true) {
        M.print(name,"gen");
    }
}
//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_msp(const Matrix<Number>& M) const {
    if (msp==true) {
        M.print(name,"msp");
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_tri() const{
    if (tri==true) {
        string file_name = name+".tri";
        ofstream out(file_name.c_str());

        const vector< pair<vector<libQnormaliz::key_t>,Number> >& Tri = Result->getTriangulation();
        typename vector< pair<vector<libQnormaliz::key_t>,Number> >::const_iterator tit = Tri.begin();        
        const vector<vector<bool> >& Dec = Result->isComputed(QConeProperty::ConeDecomposition) ?
                Result->getOpenFacets() : vector<vector<bool> >();
        typename vector< vector<bool> >::const_iterator idd = Dec.begin();

        out << Tri.size() << endl;
        size_t nr_extra_entries=1;
        if (Result->isComputed(QConeProperty::ConeDecomposition))
            nr_extra_entries+=Result->getSublattice().getRank();
        out << Result->getSublattice().getRank()+nr_extra_entries << endl; //works also for empty list

        for(; tit != Tri.end(); ++tit) {
            for (size_t i=0; i<tit->first.size(); i++) {
                out << tit->first[i] +1 << " ";
            }
            out << "   " << tit->second;
            if(Result->isComputed(QConeProperty::ConeDecomposition)){
                out << "   ";
                for (size_t i=0; i<tit->first.size(); i++) {
                    out << " " << (*idd)[i];
                }                
                idd++;
            }
            out << endl;
        }
        if (Result->isTriangulationNested()) out << "nested" << endl;
        else out << "plain" << endl;
        if (Result->isTriangulationPartial()) out << "partial" << endl;
        out.close();
    }
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_matrix_ht1(const Matrix<Number>& M) const{
    if (ht1==true) {
        M.print(name,"ht1");
    }
}

//---------------------------------------------------------------------------

string is_maximal(long a, long b) {
    return (a == b) ? " (maximal)" : "";
}

//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_inv_file() const{
    if (inv==true) {//printing .inv file
        
        string name_open=name+".inv";                              //preparing output files
        const char* file=name_open.c_str();
        ofstream inv(file);

        if (Result->isComputed(QConeProperty::VerticesOfPolyhedron)) {
            inv << "integer number_vertices_polyhedron = "
                << Result->getNrVerticesOfPolyhedron() << endl;
        }
        if (Result->isComputed(QConeProperty::ExtremeRays)) {
            size_t nr_ex_rays = Result->getNrExtremeRays();
            inv<<"integer number_extreme_rays = "<<nr_ex_rays<<endl;
        }
        if (Result->isComputed(QConeProperty::MaximalSubspace)) {
            size_t dim_max_subspace = Result->getDimMaximalSubspace();
            inv<<"integer dim_max_subspace = "<<dim_max_subspace<<endl;
        }

        inv << "integer embedding_dim = " << dim << endl;
        if (!homogeneous){
            if (Result->isComputed(QConeProperty::AffineDim))
                inv << "integer affine_dim_polyhedron = " << Result->getAffineDim() << endl;
            if (Result->isComputed(QConeProperty::RecessionRank))
                inv << "integer recession_rank = "  << Result->getRecessionRank() << endl;
        }
        if (Result->isComputed(QConeProperty::SupportHyperplanes)) { 
            inv<<"integer number_support_hyperplanes = "<<Result->getNrSupportHyperplanes()<<endl;
        }
        if (Result->isComputed(QConeProperty::TriangulationSize)) {
            inv << "integer size_triangulation = " << Result->getTriangulationSize() << endl;
        }
        if (Result->isComputed(QConeProperty::TriangulationDetSum)) {
            inv << "integer sum_dets = " << Result->getTriangulationDetSum() << endl;
        }
        
        if (!Result->isComputed(QConeProperty::Dehomogenization)) {
            inv << "boolean inhomogeneous = false" << endl;
        }
        else {
            inv << "boolean inhomogeneous = true" << endl;
            vector<Number> Linear_Form = Result->getDehomogenization();
            inv << "vector " << Linear_Form.size()
                << " dehomogenization = " << Linear_Form;
        }
        
        inv.close();
    }
}


//---------------------------------------------------------------------------

template<typename Number, typename NumberField>
void Output<Number, NumberField>::write_files() const {
    vector<libQnormaliz::key_t> rees_ideal_key;

    if (esp && Result->isComputed(QConeProperty::SupportHyperplanes) && Result->isComputed(QConeProperty::Sublattice)) {
        //write the suport hyperplanes of the full dimensional cone
        const Sublattice_Representation<Number>& BasisChange = Result->getSublattice();
        Matrix<Number> Support_Hyperplanes_Full_Cone = BasisChange.to_sublattice_dual(Result->getSupportHyperplanesMatrix());
        // Support_Hyperplanes_Full_Cone.print(name,"esp");
        string esp_string = name+".esp";
        const char* esp_file = esp_string.c_str();
        ofstream esp_out(esp_file);
        Support_Hyperplanes_Full_Cone.print(esp_out);
        esp_out << "inequalities" << endl;
        if (Result->isComputed(QConeProperty::Grading)) {
            esp_out << 1 << endl << Result->getRank() << endl;
        }
        if (Result->isComputed(QConeProperty::Dehomogenization)) {
            esp_out << 1 << endl << Result->getRank() << endl;
            esp_out << BasisChange.to_sublattice_dual(Result->getDehomogenization());
            esp_out << "dehomogenization" << endl;
        }
        esp_out.close();
    }
    if (tgn && Result->isComputed(QConeProperty::Generators))
        Result->getGeneratorsMatrix().print(name,"tgn");
    if (tri && Result->isComputed(QConeProperty::Triangulation)) {     //write triangulation
        write_tri();
    }

    if (out==true) {  //printing .out file
        string name_open=name+".out";                              //preparing output files
        const char* file=name_open.c_str();
        ofstream out(file);
        if(out.fail()){
            errorOutput() << "Cannot write to output file." << endl;
            exit(1);
        }

        // write "header" of the .out file
        
        write_renf(out);

        if (Result->isComputed(QConeProperty::VerticesOfPolyhedron)) {
            out << Result->getNrVerticesOfPolyhedron() <<" vertices of polyhedron" << endl;
        }
        if (Result->isComputed(QConeProperty::ExtremeRays)) {
            out << Result->getNrExtremeRays() <<" extreme rays" << of_cone << endl;
        }
        if (Result->isComputed(QConeProperty::SupportHyperplanes)) {
            out << Result->getNrSupportHyperplanes() <<" support hyperplanes"
                << of_polyhedron << endl;
        }
        out<<endl;

        out << "embedding dimension = " << dim << endl;
        if (homogeneous) {
            if (Result->isComputed(QConeProperty::Sublattice)) {
                auto rank = Result->getRank();
                out << "rank = "<< rank << is_maximal(rank,dim) << endl;
            }
        } else { // now inhomogeneous case
            if (Result->isComputed(QConeProperty::AffineDim))
                out << "affine dimension of the polyhedron = "
                    << Result->getAffineDim() << is_maximal(Result->getAffineDim(),dim-1) << endl;
            if (Result->isComputed(QConeProperty::RecessionRank))
                out << "rank of recession cone = "  << Result->getRecessionRank() << endl;
        }
        
        if(Result->isComputed(QConeProperty::MaximalSubspace)){
            size_t dim_max_subspace=Result->getDimMaximalSubspace();
            if(dim_max_subspace>0)
                out << "dimension of maximal subspace = " << dim_max_subspace << endl;      
        }
            
        out << endl;
        if (Result->isComputed(QConeProperty::TriangulationSize)) {
            out << "size of ";
            if (Result->isTriangulationNested()) out << "nested ";
            if (Result->isTriangulationPartial()) out << "partial ";
            out << "triangulation   = " << Result->getTriangulationSize() << endl;
        }
        if (Result->isComputed(QConeProperty::TriangulationDetSum)) {
            out << "resulting sum of |det|s = " << Result->getTriangulationDetSum() << endl;
        }
        if (Result->isComputed(QConeProperty::TriangulationSize)) {
            out << endl;
        }
        if ( Result->isComputed(QConeProperty::Dehomogenization) ) {
            out << "dehomogenization:" << endl
                << Result->getDehomogenization() << endl;
        }
        
        if ( Result->isComputed(QConeProperty::Volume)) {
            
            out << "volume (lattice normalized) = "<< Result->getVolume() << endl;
            out << "volume (Euclidean) = "<< Result->getEuclideanVolume() << endl;
            out << endl;
        }
 
        out << "***********************************************************************"
            << endl << endl;
        string module_generators_name=" lattice points in polytope";
        if (Result->isComputed(QConeProperty::ModuleGenerators)) {
            out << Result->getNrModuleGenerators() << module_generators_name <<  ":" << endl;
            Result->getModuleGeneratorsMatrix().pretty_print(out);
            out << endl;
            write_matrix_ht1(Result->getModuleGeneratorsMatrix());
        }
        
        if (Result->isComputed(QConeProperty::Deg1Elements)) {
            out << Result->getNrDeg1Elements() << module_generators_name <<  ":" << endl;
            Result->getDeg1ElementsMatrix().pretty_print(out);
            out << endl;
            write_matrix_ht1(Result->getDeg1ElementsMatrix());
        }
        
        
        if (Result->isComputed(QConeProperty::VerticesOfPolyhedron)) {
            out << Result->getNrVerticesOfPolyhedron() <<" vertices of polyhedron:" << endl;
            Result->getVerticesOfPolyhedronMatrix().pretty_print(out);
            out << endl;
        }
        if (Result->isComputed(QConeProperty::ExtremeRays)) {
            out << Result->getNrExtremeRays() << " extreme rays" << of_cone << ":" << endl;
            Result->getExtremeRaysMatrix().pretty_print(out);
            out << endl;
            if (ext) {
                // for the .gen file we append the vertices of polyhedron if there are any
                if (Result->isComputed(QConeProperty::VerticesOfPolyhedron)) {
                    Matrix<Number> Extreme_Rays(Result->getExtremeRaysMatrix());
                    Extreme_Rays.append(Result->getVerticesOfPolyhedronMatrix());
                    write_matrix_ext(Extreme_Rays);
                } else {
                    write_matrix_ext(Result->getExtremeRaysMatrix());
                }
            }
        }
        
        if(Result->isComputed(QConeProperty::MaximalSubspace) && Result->getDimMaximalSubspace()>0){
            out << Result->getDimMaximalSubspace() <<" basis elements of maximal subspace:" << endl;
            Result->getMaximalSubspaceMatrix().pretty_print(out);
            out << endl;
            if(msp)
                write_matrix_msp(Result->getMaximalSubspaceMatrix());            
        }

        //write constrains (support hyperplanes, congruences, equations)

        if (Result->isComputed(QConeProperty::SupportHyperplanes)) {
            const Matrix<Number>& Support_Hyperplanes = Result->getSupportHyperplanesMatrix();
            out << Support_Hyperplanes.nr_of_rows() <<" support hyperplanes" 
                << of_polyhedron << ":" << endl;
            Support_Hyperplanes.pretty_print(out);
            out << endl;
        }
        if (Result->isComputed(QConeProperty::Sublattice)) {
            const Sublattice_Representation<Number>& BasisChange = Result->getSublattice();
            //equations
            const Matrix<Number>& Equations = BasisChange.getEquationsMatrix();
            size_t nr_of_equ = Equations.nr_of_rows();
            if (nr_of_equ > 0) {
                out << nr_of_equ <<" equations:" <<endl;
                Equations.pretty_print(out);
                out << endl;
            }
            
            //lattice
            const Matrix<Number>& LatticeBasis = BasisChange.getEmbeddingMatrix();
            size_t nr_of_latt = LatticeBasis.nr_of_rows();
            if (nr_of_latt < dim) {
                out << nr_of_latt <<" basis elements of subspace:" <<endl;
                LatticeBasis.pretty_print(out);
                out << endl;
            }
            if(lat)
                write_matrix_lat(LatticeBasis);
            

            if (cst && Result->isComputed(QConeProperty::SupportHyperplanes)) {
                const Matrix<Number>& Support_Hyperplanes = Result->getSupportHyperplanesMatrix();
                string cst_string = name+".cst";
                const char* cst_file = cst_string.c_str();
                ofstream cst_out(cst_file);
    
                Support_Hyperplanes.print(cst_out);
                cst_out<<"inequalities"<<endl;
                Equations.print(cst_out);
                cst_out<<"equations"<<endl;

                if (Result->isComputed(QConeProperty::Dehomogenization)) {
                    cst_out << 1 << endl << dim << endl;
                    cst_out << Result->getDehomogenization();
                    cst_out << "dehomogenization" << endl;
                }
                cst_out.close();
            }
        }
        
        out.close();
    }

    write_inv_file();
}

