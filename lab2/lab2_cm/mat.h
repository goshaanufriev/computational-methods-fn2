#pragma once


double *createMat(const int &rows, const int &cols, std::ifstream &fin) {
    double *mat = new double[rows * cols];
    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j)
            fin >> mat[i * cols + j];
    }
    return mat;
}

double *idMat(const int &n) {
    double *mat = new double[n * n];
    for (auto i = 0; i < n; ++i) {
        for (auto j = 0; j < n; ++j) {
            if (i == j)
                mat[i + n * i] = 1;
            else
                mat[j + n * i] = 0;
        }
    }
    return mat;
}

double* kMat(double*& A, const int& rows, const int& cols, const double& k)
{
    double* B = new double[rows * cols];
    for (auto i = 0; i < rows * cols; ++i)
    {
        B[i] = k * A[i];
    }
    return B;
}

void multMat(double *&A, double *&B, const int &n1, const int &m1, const int &n2, const int &m2) {
    double *mat = new double[n1 * m2];
    for (auto i = 0; i < n1; ++i) {
        for (auto k = 0; k < m2; ++k) {
            mat[i * m2 + k] = 0;
            for (auto j = 0; j < n2; ++j)
                mat[i * m2 + k] += A[i * m1 + j] * B[j * m2 + k];
        }
    }
    for (auto i = 0; i < n1 * m2; ++i)
        B[i] = mat[i];
    delete[] mat;
}

double *addMat(double *&A, double *B, const int &rows, const int &cols) {
    double *C = new double[rows * cols];
    for (auto i = 0; i < rows * cols; ++i) {
        C[i] = A[i] + B[i];
    }
    return C;
}

double *subMat(double *&A, double *B, const int &rows, const int &cols) {
    double *C = new double[rows * cols];
    for (auto i = 0; i < rows * cols; ++i) {
        C[i] = A[i] - B[i];
    }
    return C;
}

double* copyMat(double* A, const int& rows, const int& cols) {
    double* B = new double[rows * cols];
    for (auto i = 0; i < rows * cols; ++i) {
        B[i] = A[i];
    }
    return B;
}

double* invDiagMat(double*& A, const int& n)
{
    double* B = new double[n];
    for (auto i = 0; i < n; ++i)
    {
        B[i] = 1 / A[i];
    }
    return B;
}

void addToDiagMat(double*& D, double*& A, const int& n)
{
    for (auto i = 0; i < n; ++i)
        A[i * n + i] += D[i];
}

double *rotMat(const int &i, const int &j, double *&A, int &n) {
    double *T = idMat(n);
    double c = A[i * n + i];
    double s = A[j * n + i];
    if (c * c + s * s != 0) {
        T[n * i + i] = c / sqrt(c * c + s * s);
        T[n * i + j] = s / sqrt(c * c + s * s);
        T[n * j + i] = -s / sqrt(c * c + s * s);
        T[n * j + j] = c / sqrt(c * c + s * s);
    }
    return T;
}

double cubeNorm(double *&A, const int &rows, const int &cols = 1) {
    double norm = 0;
    if (cols != 1) {
        double *sums = new double[rows];
        for (auto i = 0; i < rows; ++i) {
            sums[i] = 0;
            for (auto j = 0; j < cols; ++j) {
                sums[i] += std::fabs(A[i * rows + j]);
            }
            if (sums[i] > norm)
                norm = sums[i];
        }
    } else {
        for (auto i = 0; i < rows; ++i) {
            if (std::fabs(A[i]) > norm)
                norm = std::fabs(A[i]);
        }
    }
    return norm;
}

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
