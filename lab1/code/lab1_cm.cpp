#include <iostream>
#include <fstream>

#include "gauss.h"
//#include "qr.h"

double* createMat(const int& rows, const int& cols, std::ifstream& fin)
{
    double* mat = new double[rows * cols];
    for (auto i = 0; i < rows; ++i)
    {
        for (auto j = 0; j < cols; ++j)
            fin >> mat[i * cols + j];
    }
    return mat;
}

int main()
{
    std::ifstream fin("input.txt");
    int n;
    fin >> n;

    double* A = createMat(n, n, fin);
    double* b = createMat(n, 1, fin);

    double* X = gauss(n, A, b);

    for (auto i = 0; i < n; ++i)
    {
        std::cout << X[i] << "\n";
    }
    
    return 0;
}