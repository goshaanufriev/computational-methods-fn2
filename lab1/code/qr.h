#pragma once

double* qr(int& n, double*& A, double*& b)
{
    for (auto i = 0; i < n; ++i)
    {
        for (auto j = i + 1; j < n; ++j)
        {
            double* T = rotMat(i, j, A, n);
            multMat(T, A, n, n, n, n);
            multMat(T, b, n, n, n, 1);
            delete[] T;
        }
    }
    double* X = gauss(n, A, b);
    return X;
}