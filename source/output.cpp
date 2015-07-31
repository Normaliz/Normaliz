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

#include "output.h"
#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/map_operations.h"

//---------------------------------------------------------------------------

template<typename Integer>
Output<Integer>::Output(){
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
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_lattice_ideal_input(bool value){
    lattice_ideal_input=value;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::read() const{
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

template<typename Integer>
void Output<Integer>::set_name(const string& n){
    name=n;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::setCone(Cone<Integer> & C) {
    this->Result = &C;
    dim = Result->getEmbeddingDim();
    rank = Result->getRank();
    homogeneous = !Result->isInhomogeneous();
    if (homogeneous) {
        of_cone       = "";
        of_monoid     = "";
        of_polyhedron = "";
    } else {
        of_cone       = " of recession cone";
        of_monoid     = " of recession monoid";
        of_polyhedron = " of polyhedron";
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_out(const bool& flag){
    out=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_inv(const bool& flag){
    inv=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_ext(const bool& flag){
    ext=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_esp(const bool& flag){
    esp=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_typ(const bool& flag){
    typ=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_egn(const bool& flag){
    egn=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_gen(const bool& flag){
    gen=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_cst(const bool& flag){
    cst=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_tri(const bool& flag) {
    tri=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_tgn(const bool& flag) {
    tgn=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_ht1(const bool& flag) {
    ht1=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_dec(const bool& flag) {
    dec=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_mod(const bool& flag) {
    mod=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_lat(const bool& flag) {
    lat=flag;
}
//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_extra_files(){
    out=true;
    inv=true;
    gen=true;
    cst=true;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_all_files(){
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
}



//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_ext(const Matrix<Integer>& M) const{
    if (ext==true) {
        M.print(name,"ext");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_mod(const Matrix<Integer>& M) const{
    if (mod==true) {
        M.print(name,"mod");
    }
}


//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_lat(const Matrix<Integer>& M) const{
    if (ext==true) {
        M.print(name,"lat");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_esp(const Matrix<Integer>& M) const{
    if (esp==true) {
        M.print(name,"esp");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_typ(const Matrix<Integer>& M) const{
    if (typ==true) {
        M.print(name,"typ");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_egn(const Matrix<Integer>& M) const {
    if (egn==true) {
        M.print(name,"egn");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_gen(const Matrix<Integer>& M) const {
    if (gen==true) {
        M.print(name,"gen");
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_tri() const{
    if (tri==true) {
        string file_name = name+".tri";
        ofstream out(file_name.c_str());

        const vector< pair<vector<libnormaliz::key_t>,Integer> >& Tri = Result->getTriangulation();
        typename vector< pair<vector<libnormaliz::key_t>,Integer> >::const_iterator tit = Tri.begin();

        out << Tri.size() << endl;
        out << Result->getSublattice().getRank()+1 << endl; //works also for empty list

        for(; tit != Tri.end(); ++tit) {
            for (size_t i=0; i<tit->first.size(); i++) {
                out << tit->first[i]+1 << " ";
            }
            out << tit->second << endl;
        }
        out.close();
    }
}

//---------------------------------------------------------------------------


template<typename Integer>
void Output<Integer>::write_Stanley_dec() const {
    if (dec && Result->isComputed(ConeProperty::StanleyDec)) {
        ofstream out((name+".dec").c_str());

        if (Result->isComputed(ConeProperty::InclusionExclusionData)) {
            const vector< pair<vector<libnormaliz::key_t>, long> >& InExData = Result->getInclusionExclusionData();
            out << "in_ex_data" << endl;
            out << InExData.size() << endl;
            for (size_t i=0; i<InExData.size(); ++i) {
                out << InExData[i].first.size() << " ";
                for (size_t j=0; j<InExData[i].first.size(); ++j) {
                    out << InExData[i].first[j] << " ";
                }
                out << InExData[i].second << endl;  
            }
        }

        out << "Stanley_dec" << endl;
        const list<STANLEYDATA<Integer> >& StanleyDec = Result->getStanleyDec();
        typename list<STANLEYDATA<Integer> >::const_iterator S = StanleyDec.begin();
        size_t i;

        out << StanleyDec.size() << endl; 
        for (; S!=StanleyDec.end(); ++S) {
            for (i=0; i<rank; ++i)
                out << S->key[i]+1 <<" ";
            out << endl;
            S->offsets.print(out);
            out << endl;
        }
        out.close();
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_ht1(const Matrix<Integer>& M) const{
    if (ht1==true) {
        M.print(name,"ht1");
    }
}

//---------------------------------------------------------------------------

string is_maximal(long a, long b) {
    return (a == b) ? " (maximal)" : "";
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_inv_file() const{
    if (inv==true) {//printing .inv file
        size_t i;
        string name_open=name+".inv";                              //preparing output files
        const char* file=name_open.c_str();
        ofstream inv(file);

        if (Result->isComputed(ConeProperty::ModuleGenerators)) {
            inv << "integer number_module_generators = "
                << Result->getNrModuleGenerators() << endl;
        }
        if (Result->isComputed(ConeProperty::HilbertBasis)) {
            inv<<"integer hilbert_basis_elements = "<<Result->getNrHilbertBasis()<<endl;
        }

        if (Result->isComputed(ConeProperty::VerticesOfPolyhedron)) {
            inv << "integer number_vertices_polyhedron = "
                << Result->getNrVerticesOfPolyhedron() << endl;
        }
        if (Result->isComputed(ConeProperty::ExtremeRays)) {
            size_t nr_ex_rays = Result->getNrExtremeRays();
            inv<<"integer number_extreme_rays = "<<nr_ex_rays<<endl;
        }
        if (Result->isComputed(ConeProperty::ModuleGeneratorsOfIntegralClosure)) {
            inv<<"integer mod_gen_intcl = "<<Result->getNrModuleGeneratorsOfIntegralClosure()<<endl;
        }

        inv << "integer embedding_dim = " << dim << endl;
        if (homogeneous) {
            inv << "integer rank = " << rank << endl;
            inv << "integer external_index = " << Result->getSublattice().getExternalIndex() << endl;
        } else {
            if (Result->isComputed(ConeProperty::AffineDim))
                inv << "integer affine_dim_polyhedron = " << Result->getAffineDim() << endl;
            if (Result->isComputed(ConeProperty::RecessionRank))
                inv << "integer recession_rank = "  << Result->getRecessionRank() << endl;
        }
        
        if(Result->isComputed(ConeProperty::OriginalMonoidGenerators)){
            inv << "integer internal_index = " << Result->getIndex() << endl;           
        }
        if (Result->isComputed(ConeProperty::SupportHyperplanes)) { 
            inv<<"integer number_support_hyperplanes = "<<Result->getNrSupportHyperplanes()<<endl;
        }
        if (Result->isComputed(ConeProperty::TriangulationSize)) {
            inv << "integer size_triangulation = " << Result->getTriangulationSize() << endl;
        }
        if (Result->isComputed(ConeProperty::TriangulationDetSum)) {
            inv << "integer sum_dets = " << Result->getTriangulationDetSum() << endl;
        }

        if (!Result->isComputed(ConeProperty::Dehomogenization)) {
            inv << "boolean inhomogeneous = false" << endl;
        }
        else {
            inv << "boolean inhomogeneous = true" << endl;
            vector<Integer> Linear_Form = Result->getDehomogenization();
            inv << "vector " << Linear_Form.size()
                << " dehomogenization = " << Linear_Form;
        }
        if (Result->isComputed(ConeProperty::Grading)==false) {
            if (Result->isComputed(ConeProperty::ExtremeRays)) {
                inv<<"boolean graded = "<<"false"<<endl;
            }
        }
        else {
            inv<<"boolean graded = "<<"true"<<endl;
            if (Result->isComputed(ConeProperty::Deg1Elements)) {
                inv<<"integer degree_1_elements = "<<Result->getNrDeg1Elements()<<endl;
            }
            vector<Integer> Linear_Form = Result->getGrading();
            inv<<"vector "<<Linear_Form.size()<<" grading = ";
            for (i = 0; i < Linear_Form.size(); i++) {
                inv<<Linear_Form[i]<<" ";
            }
            inv<<endl;
            inv <<"integer grading_denom = " << Result->getGradingDenom() << endl;
            if (Result->isComputed(ConeProperty::Multiplicity)){
                mpq_class mult = Result->getMultiplicity();
                inv <<"integer multiplicity = "      << mult.get_num() << endl;
                inv <<"integer multiplicity_denom = "<< mult.get_den() << endl;
            }
            if (Result->isComputed(ConeProperty::HilbertSeries)) {
                const HilbertSeries& HS = Result->getHilbertSeries();
                vector<mpz_class> HSnum = HS.getNum();
                inv <<"vector "<< HSnum.size() <<" hilbert_series_num = ";
                inv << HSnum;
                vector<denom_t> HSdenom = to_vector(HS.getDenom());
                inv <<"vector "<< HSdenom.size() <<" hilbert_series_denom = ";
                inv << HSdenom;

                vector< vector <mpz_class> > hqp = HS.getHilbertQuasiPolynomial();
                if (hqp.size()>0) {
                    inv << "matrix " << hqp.size() << " " << hqp[0].size()
                        << " hilbert_quasipolynomial = ";
                    inv << endl << hqp;
                    inv << "integer hilbert_quasipolynomial_denom = " 
                        << HS.getHilbertQuasiPolynomialDenom() << endl;
                }
            }
        }
        if (Result->isComputed(ConeProperty::IsReesPrimary)) {
            if (Result->isReesPrimary()) {
                inv<<"boolean primary = true"<<endl;
            } else {
                inv<<"boolean primary = false"<<endl;
            }
        }
        if (Result->isComputed(ConeProperty::ReesPrimaryMultiplicity)) {
            inv<<"integer ideal_multiplicity = "<<Result->getReesPrimaryMultiplicity()<<endl;
        }
        
        if (Result->isComputed(ConeProperty::ClassGroup)) {
            inv <<"vector "<< Result->getClassGroup().size() <<" class_group = "<< Result->getClassGroup();
        }
        
        inv.close();
    }
}


//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_files() const {
    const Sublattice_Representation<Integer>& BasisChange = Result->getSublattice();
    size_t i, nr;
    const Matrix<Integer>& Generators = Result->getGeneratorsMatrix();
    const Matrix<Integer>& Support_Hyperplanes = Result->getSupportHyperplanesMatrix();
    vector<libnormaliz::key_t> rees_ideal_key;

    if (esp && Result->isComputed(ConeProperty::SupportHyperplanes)) {
        //write the suport hyperplanes of the full dimensional cone
        Matrix<Integer> Support_Hyperplanes_Full_Cone = BasisChange.to_sublattice_dual(Support_Hyperplanes);
        Support_Hyperplanes_Full_Cone.print(name,"esp");
    }
    if (tgn)
        Generators.print(name,"tgn");
    if (tri && Result->isComputed(ConeProperty::Triangulation)) {     //write triangulation
        write_tri();
    }

    if (out==true) {  //printing .out file
        string name_open=name+".out";                              //preparing output files
        const char* file=name_open.c_str();
        ofstream out(file);

        // write "header" of the .out file
        size_t nr_orig_gens = Result->getNrOriginalMonoidGenerators();
        if (lattice_ideal_input) {
            out << nr_orig_gens <<" original generators of the toric ring"<<endl;
        }
        if (Result->isComputed(ConeProperty::ModuleGenerators)) {
            out << Result->getNrModuleGenerators() <<" module generators" << endl;
        }
        if (Result->isComputed(ConeProperty::HilbertBasis)) {
            out << Result->getNrHilbertBasis() <<" Hilbert basis elements"
                << of_monoid << endl;
        }
        if (homogeneous && Result->isComputed(ConeProperty::Deg1Elements)) {
            out << Result->getNrDeg1Elements() <<" Hilbert basis elements of degree 1"<<endl;
        }
        if (Result->isComputed(ConeProperty::IsReesPrimary)
            && Result->isComputed(ConeProperty::HilbertBasis)) {
            const Matrix<Integer>& Hilbert_Basis = Result->getHilbertBasisMatrix();
            nr = Hilbert_Basis.nr_of_rows();
            for (i = 0; i < nr; i++) {
                if (Hilbert_Basis.read(i,dim-1)==1) {
                    rees_ideal_key.push_back(i);
                 }
            }
            out << rees_ideal_key.size() <<" generators of integral closure of the ideal"<<endl;
        }
        if (Result->isComputed(ConeProperty::VerticesOfPolyhedron)) {
            out << Result->getNrVerticesOfPolyhedron() <<" vertices of polyhedron" << endl;
        }
        if (Result->isComputed(ConeProperty::ExtremeRays)) {
            out << Result->getNrExtremeRays() <<" extreme rays" << of_cone << endl;
        }
        if(Result->isComputed(ConeProperty::ModuleGeneratorsOfIntegralClosure)) {
            out << Result->getNrModuleGeneratorsOfIntegralClosure() <<" module generators of integral closure/original monoid" << endl;    
        }
        if (Result->isComputed(ConeProperty::SupportHyperplanes)) {
            out << Result->getNrSupportHyperplanes() <<" support hyperplanes"
                << of_polyhedron << endl;
        }
        out<<endl;
        if (Result->isComputed(ConeProperty::ExcludedFaces)) {
            out << Result->getNrExcludedFaces() <<" excluded faces"<<endl;
            out << endl;
        }
        out << "embedding dimension = " << dim << endl;
        if (homogeneous) {
            out << "rank = "<< rank << is_maximal(rank,dim) << endl;
            //out << "index E:G = "<< BasisChange.get_index() << endl;
            out << "external index = "<< BasisChange.getExternalIndex() << endl;
        } else { // now inhomogeneous case
            if (Result->isComputed(ConeProperty::AffineDim))
                out << "affine dimension of the polyhedron = "
                    << Result->getAffineDim() << is_maximal(Result->getAffineDim(),dim-1) << endl;
            if (Result->isComputed(ConeProperty::RecessionRank))
                out << "rank of recession monoid = "  << Result->getRecessionRank() << endl;
        }
        
        if(Result->isComputed(ConeProperty::OriginalMonoidGenerators)){
            out << "internal index = " << Result->getIndex() << endl;           
        }
        
        if (homogeneous && Result->isComputed(ConeProperty::IsIntegrallyClosed)) {
            if (Result->isIntegrallyClosed()) {
                out << "original monoid is integrally closed"<<endl;
            } else {
                out << "original monoid is not integrally closed"<<endl;
            }
        }
        out << endl;
        if (Result->isComputed(ConeProperty::TriangulationSize)) {
            out << "size of triangulation   = " << Result->getTriangulationSize() << endl;
        }
        if (Result->isComputed(ConeProperty::TriangulationDetSum)) {
            out << "resulting sum of |det|s = " << Result->getTriangulationDetSum() << endl;
        }
        if (Result->isComputed(ConeProperty::TriangulationSize)) {
            out << endl;
        }
        if ( Result->isComputed(ConeProperty::Dehomogenization) ) {
            out << "dehomogenization:" << endl
                << Result->getDehomogenization() << endl;
        }
        if ( Result->isComputed(ConeProperty::Grading) ) {
            out << "grading:" << endl
                << Result->getGrading();
            Integer denom = Result->getGradingDenom();
            if (denom != 1) {
                out << "with denominator = " << denom << endl;
            }
            out << endl;
            if (homogeneous && Result->isComputed(ConeProperty::ExtremeRays)) {
                out << "degrees of extreme rays:"<<endl;
                map<Integer,long> deg_count;
                vector<Integer> degs = Result->getExtremeRaysMatrix().MxV(Result->getGrading());
                for (i=0; i<degs.size(); ++i) {
                    deg_count[degs[i]/denom]++;
                }
                out << deg_count;
            }
        }
        else if (Result->isComputed(ConeProperty::IsDeg1ExtremeRays)) {
            if ( !Result->isDeg1ExtremeRays() ) {
                out << "No implicit grading found" << endl;
            }
        }
        out<<endl;
        if (homogeneous && Result->isComputed(ConeProperty::IsDeg1HilbertBasis)
          && Result->isDeg1ExtremeRays() ) {
            if (Result->isDeg1HilbertBasis()) {
                out << "Hilbert basis elements are of degree 1";
            } else {
                out << "Hilbert basis elements are not of degree 1";
            }
            out<<endl<<endl;
        }
        if ( Result->isComputed(ConeProperty::ModuleRank) ) {
            out << "module rank = "<< Result->getModuleRank() << endl;
        }
        if ( Result->isComputed(ConeProperty::Multiplicity) ) {
            out << "multiplicity = "<< Result->getMultiplicity() << endl;
        }
        if ( Result->isComputed(ConeProperty::ModuleRank) || Result->isComputed(ConeProperty::Multiplicity)) {
            out << endl;
        }
        
        if ( Result->isComputed(ConeProperty::HilbertSeries) ) {
            const HilbertSeries& HS = Result->getHilbertSeries();
            out << "Hilbert series:" << endl << HS.getNum();
            map<long, long> HS_Denom = HS.getDenom();
            long nr_factors = 0;
            for (map<long, long>::iterator it = HS_Denom.begin(); it!=HS_Denom.end(); ++it) {
                nr_factors += it->second;
            }
            out << "denominator with " << nr_factors << " factors:" << endl;
            out << HS.getDenom();
            out << endl;
            if (HS.getShift() != 0) {
                out << "shift = " << HS.getShift() << endl << endl;
            }
            out << "degree of Hilbert Series as rational function = "
                << HS.getDegreeAsRationalFunction() << endl << endl;

            long period = HS.getPeriod();
            if (period == 1) {
                out << "Hilbert polynomial:" << endl;
                out << HS.getHilbertQuasiPolynomial()[0];
                out << "with common denominator = ";
                out << HS.getHilbertQuasiPolynomialDenom();
                out << endl<< endl;
            } else {
                // output cyclonomic representation
                out << "Hilbert series with cyclotomic denominator:" << endl;
                out << HS.getCyclotomicNum();
                out << "cyclotomic denominator:" << endl;
                out << HS.getCyclotomicDenom();
                out << endl;
                // Hilbert quasi-polynomial
                vector< vector<mpz_class> > hilbert_quasi_poly = HS.getHilbertQuasiPolynomial();
                if (hilbert_quasi_poly.size() > 0) { // == 0 means not computed
                    out<<"Hilbert quasi-polynomial of period " << period << ":" << endl;
                    Matrix<mpz_class> HQP(hilbert_quasi_poly);
                    HQP.pretty_print(out,true);
                    out<<"with common denominator = "<<HS.getHilbertQuasiPolynomialDenom();
                }
                out << endl << endl;
            }

        }

        if (Result->isComputed(ConeProperty::IsReesPrimary)) {
            if (Result->isReesPrimary()) {
                out<<"ideal is primary to the ideal generated by the indeterminates"<<endl;
            } else {
                out<<"ideal is not primary to the ideal generated by the indeterminates"<<endl;
            }
            if (Result->isComputed(ConeProperty::ReesPrimaryMultiplicity)) {
                out<<"multiplicity of the ideal = "<<Result->getReesPrimaryMultiplicity()<<endl;
            }
            out << endl;
        }
        
        if(Result->isComputed(ConeProperty::ClassGroup)) {
            vector<Integer> ClassGroup=Result->getClassGroup();
            out << "rank of class group = " << ClassGroup[0] << endl;
            if(ClassGroup.size()==1)
                out << "class group is free" << endl << endl;
            else{
                ClassGroup.erase(ClassGroup.begin());
                out << "finite cyclic summands:" << endl;
                out << count_in_map<Integer,size_t>(ClassGroup);
                out << endl;
            }                   
        }

        out << "***********************************************************************"
            << endl << endl;


        if (lattice_ideal_input) {
            out << nr_orig_gens <<" original generators:"<<endl;
            Result->getOriginalMonoidGeneratorsMatrix().pretty_print(out);
            out << endl;
        }
        if (Result->isComputed(ConeProperty::ModuleGenerators)) {
            out << Result->getNrModuleGenerators() <<" module generators:" << endl;
            Result->getModuleGeneratorsMatrix().pretty_print(out);
            out << endl;
        }
        
        if ( Result->isComputed(ConeProperty::Deg1Elements) ) {
            const Matrix<Integer>& Hom = Result->getDeg1ElementsMatrix();
            write_matrix_ht1(Hom);
            nr=Hom.nr_of_rows();
            out<<nr<<" Hilbert basis elements of degree 1:"<<endl;
            Hom.pretty_print(out);
            out << endl;
        }
        
        if (Result->isComputed(ConeProperty::HilbertBasis)) {

            const Matrix<Integer>& Hilbert_Basis = Result->getHilbertBasisMatrix();
            
            if(!Result->isComputed(ConeProperty::Deg1Elements)){
                nr=Hilbert_Basis.nr_of_rows();
                out << nr << " Hilbert basis elements" << of_monoid << ":" << endl;
                Hilbert_Basis.pretty_print(out);
                out << endl;
            }
            else {
                nr=Hilbert_Basis.nr_of_rows()-Result->getNrDeg1Elements();
                out << nr << " further Hilbert basis elements" << of_monoid << " of higher degree:" << endl;
                Matrix<Integer> HighDeg(nr,dim);
                for(size_t i=0;i<nr;++i)
                    HighDeg[i]=Hilbert_Basis[i+Result->getNrDeg1Elements()];                               
                HighDeg.pretty_print(out);
                out << endl;                       
            }
            Matrix<Integer> complete_Hilbert_Basis(0,dim);
            if (gen || egn || typ) {
                // for these files we append the module generators if there are any
                if (Result->isComputed(ConeProperty::ModuleGenerators)) {
                    complete_Hilbert_Basis.append(Hilbert_Basis);
                    complete_Hilbert_Basis.append(Result->getModuleGeneratorsMatrix());
                    write_matrix_gen(complete_Hilbert_Basis);
                } else {
                    write_matrix_gen(Hilbert_Basis);
                }
            }
            if (egn || typ) {
                Matrix<Integer> Hilbert_Basis_Full_Cone = BasisChange.to_sublattice(Hilbert_Basis);
                if (egn) {
                    string egn_string = name+".egn";
                    const char* egn_file = egn_string.c_str();
                    ofstream egn_out(egn_file);
        
                    Hilbert_Basis_Full_Cone.print(egn_out);
                    egn_out<<"integral_closure"<<endl;
                    if (Result->isComputed(ConeProperty::Grading)) {
                        egn_out << 1 << endl << rank << endl;
                        egn_out << BasisChange.to_sublattice_dual(Result->getGrading());
                        egn_out << "grading" << endl;
                    }
                    if (Result->isComputed(ConeProperty::Dehomogenization)) {
                        egn_out << 1 << endl << rank << endl;
                        egn_out << BasisChange.to_sublattice_dual(Result->getDehomogenization());
                        egn_out << "dehomogenization" << endl;
                    }                    
                    egn_out.close();
                }    

                if (typ && homogeneous) {
                    write_matrix_typ(Hilbert_Basis_Full_Cone.multiplication(BasisChange.to_sublattice_dual(Support_Hyperplanes).transpose()));
                }
            }

            if (Result->isComputed(ConeProperty::IsReesPrimary)) {
                out << rees_ideal_key.size() <<" generators of integral closure of the ideal:"<<endl;
                Matrix<Integer> Ideal_Gens = Hilbert_Basis.submatrix(rees_ideal_key);
                Ideal_Gens.resize_columns(dim-1);
                Ideal_Gens.pretty_print(out);
                out << endl;
            }
        }
        if (Result->isComputed(ConeProperty::VerticesOfPolyhedron)) {
            out << Result->getNrVerticesOfPolyhedron() <<" vertices of polyhedron:" << endl;
            Result->getVerticesOfPolyhedronMatrix().pretty_print(out);
            out << endl;
        }
        if (Result->isComputed(ConeProperty::ExtremeRays)) {
            out << Result->getNrExtremeRays() << " extreme rays" << of_cone << ":" << endl;
            Result->getExtremeRaysMatrix().pretty_print(out);
            out << endl;
            if (ext) {
                // for the .gen file we append the vertices of polyhedron if there are any
                if (Result->isComputed(ConeProperty::VerticesOfPolyhedron)) {
                    Matrix<Integer> Extreme_Rays(Result->getExtremeRaysMatrix());
                    Extreme_Rays.append(Result->getVerticesOfPolyhedronMatrix());
                    write_matrix_ext(Extreme_Rays);
                } else {
                    write_matrix_ext(Result->getExtremeRaysMatrix());
                }
            }
        }
        
        if(Result->isComputed(ConeProperty::ModuleGeneratorsOfIntegralClosure)) {
            out << Result->getNrModuleGeneratorsOfIntegralClosure() <<" module generators of integral closure/original monoid:" << endl;
            Result->getModuleGeneratorsOfIntegralClosureMatrix().pretty_print(out);
            out << endl;
            if(mod)
                write_matrix_mod(Result->getModuleGeneratorsOfIntegralClosureMatrix());
        }

        //write constrains (support hyperplanes, congruences, equations)

        if (Result->isComputed(ConeProperty::SupportHyperplanes)) {
            out << Support_Hyperplanes.nr_of_rows() <<" support hyperplanes" 
                << of_polyhedron << ":" << endl;
            Support_Hyperplanes.pretty_print(out);
            out << endl;
        }
        if (Result->isComputed(ConeProperty::ExtremeRays)) {
            //equations
            const Matrix<Integer>& Equations = BasisChange.getEquationsMatrix();
            size_t nr_of_equ = Equations.nr_of_rows();
            if (nr_of_equ > 0) {
                out << nr_of_equ <<" equations:" <<endl;
                Equations.pretty_print(out);
                out << endl;
            }

            //congruences
            const Matrix<Integer>& Congruences = BasisChange.getCongruencesMatrix();
            size_t nr_of_cong = Congruences.nr_of_rows();
            if (nr_of_cong > 0) {
                out << nr_of_cong <<" congruences:" <<endl;
                Congruences.pretty_print(out);
                out << endl;
            }
            
            //lattice
            const Matrix<Integer>& LatticeBasis = BasisChange.getEmbeddingMatrix();
            size_t nr_of_latt = LatticeBasis.nr_of_rows();
            if (nr_of_latt < dim ||  BasisChange.getExternalIndex()!=1) {
                out << nr_of_latt <<" basis elements of lattice:" <<endl;
                LatticeBasis.pretty_print(out);
                out << endl;
            }
            if(lat)
                write_matrix_lat(LatticeBasis);
            

            //excluded faces
            const Matrix<Integer>& ExFaces = Result->getExcludedFacesMatrix();
            size_t nr_of_exfaces = 0;
            if (Result->isComputed(ConeProperty::ExcludedFaces)) {
                nr_of_exfaces = ExFaces.nr_of_rows();
                out << nr_of_exfaces <<" excluded faces:" <<endl;
                ExFaces.pretty_print(out);
                out << endl;
            }

            if(cst) {
                string cst_string = name+".cst";
                const char* cst_file = cst_string.c_str();
                ofstream cst_out(cst_file);
    
                Support_Hyperplanes.print(cst_out);
                cst_out<<"inequalities"<<endl;
                Equations.print(cst_out);
                cst_out<<"equations"<<endl;
                Congruences.print(cst_out);
                cst_out<<"congruences"<<endl;
                if (Result->isComputed(ConeProperty::ExcludedFaces)) {
                    ExFaces.print(cst_out);
                    cst_out<<"excluded_faces"<<endl;
                }
                if (Result->isComputed(ConeProperty::Grading)) {
                    cst_out << 1 << endl << dim << endl;
                    cst_out << Result->getGrading();
                    cst_out << "grading" << endl;
                }
                if (Result->isComputed(ConeProperty::Dehomogenization)) {
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
    write_Stanley_dec();
}

