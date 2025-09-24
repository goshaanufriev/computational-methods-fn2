#include <iostream>
#include <fstream>
int err{ 0 };


#include "mat.h"
#include "gauss.h"
#include "qr.h"

double* invMat(const int& n, double*& A)
{
    double* invA = new double[n * n];
    double* idA = idMat(n);
    for (auto j = 0; j < n; ++j)
    {
        double* copyA = new double[n * n];
        for (auto i = 0; i < n * n; ++i)
        {
            copyA[i] = A[i];
        }
        double* ej = new double[n];
        for (auto i = 0; i < n; ++i)
        {
            ej[i] = idA[i * n + j];
        }
        double* X = gauss(n, copyA, ej);
        for (auto i = 0; i < n; ++i)
        {
            invA[i * n + j] = X[i];
        }
        delete[] X, copyA, ej;
    }
    return invA;
}

double cubeNorm(const double*& A, const int& n)
{
    double* sums = new double[n];
    double norm = 0;
    for (auto i = 0; i < n; ++i)
    {
        sums[i] = 0;
        for (auto j = 0; j < n; ++j)
        {
            sums[i] += std::fabs(A[i * n + j]);
        }
        if (sums[i] > norm)
            norm = sums[i];
    }
    return norm;
}

double octNorm(const double*& A, const int& n)
{
    double* sums = new double[n];
    double norm = 0;
    for (auto j = 0; j < n; ++j)
    {
        sums[j] = 0;
        for (auto i = 0; i < n; ++i)
        {
            sums[j] += std::fabs(A[i * n + j]);
        }
        if (sums[j] > norm)
            norm = sums[j];
    }
    return norm;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    std::ifstream fin("input.txt");
    int n;
    fin >> n;

    double* A = createMat(n, n, fin);
    //double* b = createMat(n, 1, fin);

    fin.close();

    //double* X = qr(n, A, b);

    double* invA = invMat(n, A);
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < n; ++j)
        {
            std::cout << invA[i * n + j] << " ";
        }
        std::cout << "\n";
    }
    multMat(A, invA, n, n, n, n);
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < n; ++j)
        {
            std::cout << invA[i * n + j] << " ";
        }
        std::cout << "\n";
    }

    std::ofstream fout("output4_mod.txt");

    /*if (err == 0)
    {
        for (auto i = 0; i < n; ++i)
        {
            fout << X[i] << "\n";
        }
    }*/

    fout.close();
    
    
    return 0;
}