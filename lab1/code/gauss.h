#pragma once


float* gauss(const int& n, float*& A, float*& b)
{
    // Прямой ход метода Гаусса

    for (auto j = 0; j < n; ++j)
    {
        auto maxRow = j;
        for (auto i = j + 1; i < n; ++i)
        {
            if (std::fabs(A[n * i + j]) > std::fabs(A[n * maxRow + j]))
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
        if (std::fabs(A[j * n + j]) < 1E-10)
        {
            err = 1;
            std::cout << "Ошибка! Система вырождена!\n";
            break;
        }
        float coef = A[j + n * j];
        b[j] /= coef;
        for (auto i = j; i < n; ++i)
        {
            A[i + n * j] /= coef;
        }
        for (int i = j + 1; i < n; i++) {
            float coef = A[i * n + j];
            for (int k = j; k < n; k++) {
                A[i * n + k] -= coef * A[j * n + k];
            }
            b[i] -= coef * b[j];
        }
    }

    // Обратный ход метода Гаусса

    float* X = new float[n];
    if (err == 0)
    {
        X[n - 1] = b[n - 1];
        for (int i = n - 2; i >= 0; --i)
        {
            float s = 0;
            for (auto j = i + 1; j < n; ++j)
            {
                s += A[j + n * i] * X[j];
            }
            X[i] = b[i] - s;
        }
    }
    return X;
}
