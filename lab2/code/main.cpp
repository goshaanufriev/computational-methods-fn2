#include <iostream>
#include <fstream>
#include <cmath>

#include "mat.h"

int main() {
    std::ifstream fin("input.txt");
    int n;
    fin >> n;

    double *A = createMat(n, n, fin);
    double *b = createMat(n, 1, fin);

    // double *L = new double[n * n];
    // double *U = new double[n * n];
    // double *D = new double[n];
    //
    // for (auto i = 0; i < n; ++i) {
    //     for (auto j = 0; j < n; ++j) {
    //         if (i == j) {
    //             L[i * n + j] = 0;
    //             U[i * n + j] = 0;
    //             D[i] = A[i * n + j];
    //         } else if (i > j) {
    //             L[i * n + j] = A[i * n + j];
    //             U[i * n + j] = 0;
    //         } else {
    //             L[i * n + j] = 0;
    //             U[i * n + j] = A[i * n + j];
    //         }
    //     }
    // }

    double *C = new double[n * n];
    double *y = new double[n];

    for (auto i = 0; i < n; ++i) {
        for (auto j = 0; j < n; ++j) {
            if (i == j) {
                C[i * n + j] = 0;
                y[i] = b[i] / A[i * n + j];
            }
            else {
                C[i * n + j] = -A[i * n + j] / A[i * n + i];
            }
        }
    }

    double* x0 = new double[n];
    for (auto i = 0; i < n; ++i) {
        x0[i] = 0;
    }

    double normC = cubeNorm(C,n,n);
    std::cout << normC << "\n";
    double eps = 1e-8;
    double* cx = copyMat(x0,n,n);
    multMat(C, cx,n,n,n,1);
    double* x = addMat(cx,y,n,1);
    double* err = subMat(x,x0,n,1);
    double errNorm = cubeNorm(err,n,1);
    while (errNorm > (1 - normC) * eps / normC) {
        x0 = copyMat(x,n,1);
        cx = copyMat(x0,n,1);
        multMat(C, cx,n,n,n,1);
        x = addMat(cx,y,n,1);
        err = subMat(x,x0,n,1);
        errNorm = cubeNorm(err,n,1);
    }
    for (auto i = 0; i < n; ++i) {
        std::cout << x[i] << "\n";
    }

    return 0;
}
