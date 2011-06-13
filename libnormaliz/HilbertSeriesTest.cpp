#include "HilbertSeries.cpp"
#include <iostream>
using namespace std;
using namespace libnormaliz;

/*
template <typename Class>
ostream& operator<< (ostream& out, const vector<Class>& vec) {
    for (size_t i=0; i<vec.size(); ++i) {
		out << " " << vec[i];
	}
	out << endl;
	return out;
}
*/

int main() {
	
	vector<long> anom(3);
	anom[1]=1; anom[2]=1;
	vector<long> adenom(4);
	adenom[2]=1; adenom[3]=1;
	
	HilbertSeries<long> A(anom,adenom);

	vector<long> bnom(5);
	bnom[0]=1; bnom[2]=1; bnom[4]=1;
	vector<long> bdenom(4);
	bdenom[3]=2;
	
	HilbertSeries<long> B(bnom,bdenom);

	HilbertSeries<long> ABA = B;
	cout << "A: " << A;
	cout << "B: " << B;
	cout << "B: " << ABA;
	ABA += A;
	cout << "B+A: " << ABA;
	ABA += A;
	cout << "B+A+A: " << ABA;



	return 0;
}
