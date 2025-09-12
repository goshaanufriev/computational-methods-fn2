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

    // Прямой ход метода Гаусса

    for (auto j = 0; j < n; ++j)
    {
        /*if (j > 0)
        {
            for (auto k = j - 1, i = 0; k < n; ++k, ++i)
            {
                A[k + n * j] -= A[k + n * (j - 1)] * A[j - 1 + n * j];
                b[j + i] -= A[j - 1 + n * j] * b[j - 1];
            }
        }*/
        auto maxRow = j;
        for (auto i = j + 1; i < n; ++i)
        {
            if (std::fabs(A[j + n * i]) > std::fabs(A[j + n * maxRow]))
            {
                maxRow = i;
            }
        }
        if (maxRow != j) {
            for (int k = j; k < n; k++) {
                std::swap(A[j * n + k], A[maxRow * n + k]);
            }
            std::swap(b[j], b[maxRow]);
        }
        b[j] /= A[j + n * j];
        for (auto i = j; i < n; ++i)
        {
            A[i + n * j] /= A[j + n * j];
        }
        for (int i = j + 1; i < n; i++) {
            for (int k = j; k < n; k++) {
                A[i * n + k] -= A[i * n + j] * A[j * n + k];
            }
            b[i] -= A[i * n + j] * b[j];
        }
    }

    // Обратный ход метода Гаусса

    double* X = new double[n];
    X[n-1] = b[n-1];
    for (int i = n - 2; i >= 0; --i)
    {
        double s = 0;
        for (auto j = i + 1; j < n; ++j)
        {
            s += A[j + n * i]*X[j];
        }
        X[i] = b[i] - s;
    }

    for (auto i = 0; i < n; ++i)
    {
        std::cout << X[i] << "\n";
    }
    
    return 0;
}