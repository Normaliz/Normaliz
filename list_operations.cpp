/*
 * Normaliz 2.2
 * Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
 * With contributions by Christof Soeger
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


//---------------------------------------------------------------------------

#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "simplex.h"
#include "list_operations.h"

//---------------------------------------------------------------------------

int l_read(const list< vector<Integer> >& l){
	list< vector<Integer> >::const_iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		v_read(*i);
	}
	return l.size();
}

//---------------------------------------------------------------------------

int l_read(const list< vector<int> >& l){
	list< vector<int> >::const_iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		v_read(*i);
	}
	return l.size();
}

//---------------------------------------------------------------------------

int l_read(const list< Simplex >& l){
	list< Simplex >::const_iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		(*i).read_k();
	}
	return l.size();
}

//---------------------------------------------------------------------------

int l_read(const list< int >& l){
	list< int >::const_iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		cout<<(*i)<<" ";
	}
	cout<<endl;
	return l.size();
}

//---------------------------------------------------------------------------

vector<Integer> l_multiplication(const list< vector<Integer> >& l,const vector<Integer>& v){
	register int s=l.size();
	vector<Integer> p(s);
	list< vector<Integer> >::const_iterator i;
	s=0;
	for (i =l.begin(); i != l.end(); i++) {
		p[s]=v_scalar_product(*i,v);             //maybe we loose time here?
		s++;
	}
	return p;
}

//---------------------------------------------------------------------------

list< vector<Integer> > l_list_x_matrix(const list< vector<Integer> >& l,const Matrix& M){
	list< vector<Integer> > result;
	vector<Integer> p;
	list< vector<Integer> >::const_iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		p=M.VxM(*i);
		result.push_back(p);
	}
	return result;
}
//---------------------------------------------------------------------------

void  l_cut(list<  vector<Integer> >& l, int size){
	list< vector<Integer> >::iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		(*i).resize(size);
	}
}

//---------------------------------------------------------------------------


void  l_cut_front(list<  vector<Integer> >& l, int size){
	int s,k;
	list< vector<Integer> >::iterator i;
	for (i =l.begin(); i != l.end(); i++) {
		s=(*i).size()-size;
		for (k = 0; k < size; k++) {
			(*i)[k]=(*i)[s+k];
			//			swap((*i)[k],(*i)[s+k]);  //TODO lieber nicht swap???  einfach (*i)[k]=(*i)[s+k]
		}
		(*i).resize(size);
	}
}

//---------------------------------------------------------------------------
