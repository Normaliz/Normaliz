#include "HilbertSeries.cpp"
#include <iostream>
using namespace std;
using namespace libnormaliz;

int main() {

	vector<long long> anom(3);
	anom[1]=1; anom[2]=1;
	vector<long long> adenom(4);
	adenom[2]=1; adenom[3]=1;
	
	HilbertSeries A(anom,adenom);

	vector<long long> bnom(5);
	bnom[0]=1; bnom[2]=1; bnom[4]=1;
	vector<long long> bdenom(4);
	bdenom[3]=2;
	
	HilbertSeries B(bnom,bdenom);

	HilbertSeries ABA = HilbertSeries();
	ABA += B;
	cout << "A: " << A;
	cout << "B: " << B;
	cout << "B: " << ABA;
	ABA += A;
	cout << "B+A: " << ABA;
	ABA += A;
	cout << "B+A+A: " << ABA << endl;

	ABA.simplify();
	cout << "Simpl: " << ABA;

	return 0;
}
