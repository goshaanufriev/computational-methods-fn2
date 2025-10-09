#include <iostream>
#include <fstream>
int err{ 0 };


#include "mat.h"
#include "gauss.h"
#include "qr.h"

float* invMat(const int& n, float*& A)
{
    float* invA = new float[n * n];
    float* idA = idMat(n);
    for (auto j = 0; j < n; ++j)
    {
        float* copyA = new float[n * n];
        for (auto i = 0; i < n * n; ++i)
        {
            copyA[i] = A[i];
        }
        float* ej = new float[n];
        for (auto i = 0; i < n; ++i)
        {
            ej[i] = idA[i * n + j];
        }
        float* X = gauss(n, copyA, ej);
        for (auto i = 0; i < n; ++i)
        {
            invA[i * n + j] = X[i];
        }
        delete[] X, copyA, ej;
    }
    return invA;
}

float cubeNorm(float*& A, const int& n)
{
    float* sums = new float[n];
    float norm = 0;
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

float octNorm(float*& A, const int& n)
{
    float* sums = new float[n];
    float norm = 0;
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

void print(std::ostream& out, float*& A, const int& n, const int& m = 1)
{
    for (auto i = 0; i < n; ++i)
    {
        if (m == 1)
        {
            out << A[i];
        }
        else
        {
            for (auto j = 0; j < m; ++j)
            {
                out << A[i * n + j] << " ";
            }
        }
        out << "\n";
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");

    std::ifstream fin("input_good.txt");
    int n;
    fin >> n;

    float* A = createMat(n, n, fin);
    //float* b = createMat(n, 1, fin);


    fin.close();

    float* invA = invMat(n, A);
    std::cout << "\n";
    float nA, nInvA;
    nA = octNorm(A, n);
    nInvA = octNorm(invA, n);


    std::ofstream fout("output_good_float.txt");
    fout << "Число обусловленности в октаэдрической норме: " << nA * nInvA << "\n";
    nA = cubeNorm(A, n);
    nInvA = cubeNorm(invA, n);
    fout << "Число обусловленности в кубической норме: " << nA * nInvA << "\n";
    fout << "\nA^-1*A:\n";
    multMat(invA, A, n, n, n, n);
    print(fout, A, n, n);


    /*float* X = gauss(n, A, b);
    if(err == 0)
        print(fout, X, n);*/


    fout.close();
    
    
    return 0;
}