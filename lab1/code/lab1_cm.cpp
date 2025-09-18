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

double* idMat(const int& n)
{
    double* mat = new double[n];
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = 0; j < n; ++j)
        {
            if (i == j)
                mat[i + n * i] = 1;
            else
                mat[j + n * i] = 0;
        }
    }
    return mat;
}

double*& multMat(const double*& A, double*& B, const int& n1, const int& m1, const int& n2, const int& m2)
{
    double* mat = new double[n1 * m2];
    for (auto i = 0; i < n1; ++i)
    {
        for (auto k = 0; k < m2; ++k)
        {
            mat[i * m2 + k] = 0;
            for (auto j = 0; j < n2; ++j)
                mat[i * m2 + k] += A[i * m1 + j] * B[j * m2 + k];
        }
    }
    delete[] B;
    B = mat;
    return B;
}

double* rotMat(int& i, int& j, const double*& A, int& n)
{
    double* T = idMat(n);
    double c = A[i * n + i];
    double s = A[i * n + j];
    T[n * i + i] = c / sqrt(c * c + s * s);
    T[n * i + j] = s / sqrt(c * c + s * s);
    T[n * j + i] = -s / sqrt(c * c + s * s);
    T[n * j + j] = c / sqrt(c * c + s * s);
    return T;

}

int main()
{
    setlocale(LC_ALL, "Russian");

    std::ifstream fin("input.txt");
    int n;
    fin >> n;

    double* A = createMat(n, n, fin);
    double* b = createMat(n, 1, fin);

    //double* C = new double[n * n];
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = i + 1; j < n; ++j)
        {
            double* T = rotMat(i, j, A, n);
            A = multMat(T, A, n, n, n, n);
            b = multMat(T, b, n, n, n, 1);
        }
    }
    double* X = gauss(n, A, b);

    for (auto i = 0; i < n; ++i)
    {
        std::cout << X[i] << "\n";
    }
    
    return 0;
}