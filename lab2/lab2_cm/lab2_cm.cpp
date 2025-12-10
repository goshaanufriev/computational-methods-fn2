#include <iostream>
#include <fstream>
#include <cmath>

int err{ 0 };

#include "gauss.h"
#include "mat.h"
#include "iters.h"



int main() {
	std::ifstream fin("inputVar4.txt");
	int n;
	fin >> n;

	double* A = createMat(n, n, fin);
	double* b = createMat(n, 1, fin);

	double w = 0.5;

	double* x = relax(A,b,n, w);
	
	std::ofstream fout("output.txt");
	for (auto i = 0; i < n; ++i) {
		fout << x[i] << "\n";
	}

	return 0;
}
