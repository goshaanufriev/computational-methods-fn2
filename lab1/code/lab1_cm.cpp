#include <iostream>
#include <fstream>

double* createMat(const int& rows, const int& cols, std::ifstream& fin)
{
    double* mat = new double[rows * cols];
    for (auto i = 0; i < rows; ++i)
    {
        for (auto j = 0; j < cols; ++j)
            fin >> mat[j + i * cols];
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
    
    // double* X = new double[n];

    for (auto j = 0; j < n; ++j)
    {
        auto maxInd = j + n * j;
        auto maxEl = A[maxInd];
        auto maxRow = j;
        for (auto i = j + 1; i < n; ++i)
        {
            if (abs(A[j + n * i]) > maxEl)
            {
                maxEl = A[j + n * i];
                maxInd = j + n * i;
                maxRow = i;
            }
        }
        if (maxInd != j + n * j)
        {
            for (auto i = j, k = 0; i < n; ++j, ++k)
            {
                std::swap(A[j + n * i], A[maxInd + k]);
            }
            std::swap(b[j],b[maxRow]);
        }
        for (auto i = j; i < n; ++i)
        {
            A[j + n * i] /= A[j + n * j];
        }
    }
    
    return 0;
}