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

#ifndef CONE_HELPER_H_
#define CONE_HELPER_H_

#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"

namespace libnormaliz {
using std::vector;

// determines the maximal subsets in a vector of subsets given by their indicator vectors
// result returned in is_max_subset -- must be initialized outside
// only set to false in this routine
// if a set occurs more than once, only the last instance is recognized as maximal
void maximal_subsets(const vector<vector<bool> >& ind, vector<bool>& is_max_subset);

// the project-and-lift algorithm for lattice points on a polytope

template<typename IntegerPL, typename IntegerRet> 
class ProjectAndLift {
    
    template<typename,typename> friend class ProjectAndLift;
    
    vector<Matrix<IntegerPL> > AllSupps;
    vector<vector<size_t> > AllOrders;
    
    vector<boost::dynamic_bitset<> > StartInd;
    size_t StartRank;
    
    Matrix<IntegerRet> Deg1Points;
    vector<IntegerRet> excluded_point;
    IntegerRet GD;
    
    size_t EmbDim;
    bool verbose;
    
    vector<size_t> order_supps(const Matrix<IntegerPL>& Supps);   
    bool fiber_interval(IntegerRet& MinInterval, IntegerRet& MaxInterval,
                        const vector<IntegerRet>& base_point);    
    
    void lift_point_recursively(vector<IntegerRet>& final_latt_point, 
                                const vector<IntegerRet>& latt_point_proj);    
    void lift_points_to_this_dim(Matrix<IntegerRet>& Deg1Points, const Matrix<IntegerRet>& Deg1Proj); 
    
    void find_single_point();
    void lift_points_by_generation();
    
    void compute_projections(size_t dim, vector< boost::dynamic_bitset<> >& Ind, size_t rank);
    
    void initialize(const Matrix<IntegerPL>& Supps,size_t rank);
        
    public:
        
    ProjectAndLift(const Matrix<IntegerPL>& Supps,const vector<boost::dynamic_bitset<> >& Ind,size_t rank);
    
    void set_excluded_point(const vector<IntegerRet>& excl_point);
    void set_grading_denom(const IntegerRet GradingDenom);
    void set_verbose(bool on_off);
    
    void compute(bool do_all_points=true);    
    void put_eg1Points_into(Matrix<IntegerRet>& LattPoints);    
};
    
    
template<typename IntegerPL, typename IntegerRet>
void project_and_lift_inner(Matrix<IntegerRet>& Deg1, const Matrix<IntegerPL>& Supps, 
                            vector<boost::dynamic_bitset<> >& Ind, const IntegerRet& GD, size_t rank,
                            bool verbose, bool all_points, const vector<IntegerRet>& excluded_point);

// computes c1*v1-c2*v2
template<typename Integer>
vector<Integer> FM_comb(Integer c1, const vector<Integer>& v1,Integer c2, const vector<Integer>& v2, bool& is_zero);

} //end namespace libnormaliz

#endif /* CONE_HELPER_H_ */
