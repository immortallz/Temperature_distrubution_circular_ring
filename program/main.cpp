#include "r.h"

int main()
{
	vector<double> error = 
		temperature_ring(
			100, // N ~ radial count of nodes
			360, // M ~ polar count of nodes
			0.1, // A
			0.1, // B
			1.0, // R1
			2.0, // R2
			5.0, // tau
			"analytic.txt",
			"numerical.txt",
			"error.txt");
	cout << "\nmain.cpp output:" << endl;
	for(double elem : error)
		cout << elem << endl;
	return 0;
}