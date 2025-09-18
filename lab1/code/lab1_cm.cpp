#include <iostream>
#include <fstream>
int err{ 0 };


#include "mat.h"
#include "gauss.h"
#include "qr.h"

int main()
{
    setlocale(LC_ALL, "Russian");

    std::ifstream fin("input.txt");
    int n;
    fin >> n;

    double* A = createMat(n, n, fin);
    double* b = createMat(n, 1, fin);

    double* X = gauss(n, A, b);

    if (err == 0)
    {
        for (auto i = 0; i < n; ++i)
        {
            std::cout << X[i] << "\n";
        }
    }
    
    
    return 0;
}