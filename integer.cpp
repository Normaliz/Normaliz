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

/*#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>*/
#include <algorithm>
//using namespace std;

//---------------------------------------------------------------------------

#include "integer.h"

//---------------------------------------------------------------------------

/* inline now
Integer Iabs(const Integer& a){
	if (a>=0) {
		return a;
	}
	else
		return -a;
}
*/

//---------------------------------------------------------------------------

Integer gcd(const Integer& a, const Integer& b){
	if (a==0) {
		return Iabs(b);
	}
	if (b==0) {
		return Iabs(a);
	}
	Integer q0,q1,r;
	q0=Iabs(a);
	r=Iabs(b);
	do {
		q1=r;
		r=q0%q1;
		q0=q1;
	} while (r!=0);
	return q1;
}

//---------------------------------------------------------------------------

Integer lcm(const Integer& a, const Integer& b){
	if ((a==0)||(b==0)) {
		return 0;
	}
	else
		return Iabs(a*b/gcd(a,b));
}

//---------------------------------------------------------------------------

int decimal_length(Integer a){
	int l=1;
	if (a<0) {
		a=-a;
		l++;
	}
	while((a/=10)!=0)
		l++;
	return l;
}

//---------------------------------------------------------------------------

Integer permutations(const int& a, const int& b){
	int i;
	Integer P=1;
	for (i = a+1; i <=b; i++) {
		P*=i;
	}
	return P;
}

//---------------------------------------------------------------------------

