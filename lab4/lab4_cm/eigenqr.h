#ifndef UNTITLED3_EIGENQR_H
#define UNTITLED3_EIGENQR_H

void GivensRotation(MatrixNN& A, size_t i, size_t j, double c, double s) {
    size_t n = A.getSize();

    for (size_t row = 0; row < n; row++) {
        double temp_i = c * A[row][i] + s * A[row][j];
        double temp_j = -s * A[row][i] + c * A[row][j];
        A[row][i] = temp_i;
        A[row][j] = temp_j;
    }

    for (size_t col = 0; col < n; col++) {
        double temp_i = c * A[i][col] + s * A[j][col];
        double temp_j = -s * A[i][col] + c * A[j][col];
        A[i][col] = temp_i;
        A[j][col] = temp_j;
    }
}

void toHessenberg(MatrixNN& A) {
    size_t n = A.getSize();

    for (size_t col = 0; col < n - 2; col++) {
        for (size_t row = col + 2; row < n; row++) {
            if (std::fabs(A[row][col]) > 1e-10) {
                double a = A[col + 1][col];
                double b = A[row][col];
                double r = std::sqrt(a * a + b * b);

                if (r < 1e-10) continue;

                double c = a / r;
                double s = b / r;

                GivensRotation(A, col + 1, row, c, s);
            }
        }
    }
}

void qrDecompositionHessenberg(const MatrixNN& H, MatrixNN& Q, MatrixNN& R) {
    size_t n = H.getSize();

    for (size_t i = 0; i < n - 1; i++) {
        double a = R[i][i];
        double b = R[i + 1][i];
        double r = std::sqrt(a * a + b * b);

        if (r < 1e-12) continue;

        double c = a / r;
        double s = b / r;

        for (size_t j = i; j < n; j++) {
            double temp1 = c * R[i][j] + s * R[i + 1][j];
            double temp2 = -s * R[i][j] + c * R[i + 1][j];
            R[i][j] = temp1;
            R[i + 1][j] = temp2;
        }

        for (size_t row = 0; row < n; row++) {
            double temp1 = c * Q[row][i] + s * Q[row][i + 1];
            double temp2 = -s * Q[row][i] + c * Q[row][i + 1];
            Q[row][i] = temp1;
            Q[row][i + 1] = temp2;
        }

    }
//    std::cout << "Q:" << std::endl;
//    Q.print();
//    std::cout << "R:" << std::endl;
//    R.print();
}

std::vector<double> qrAlgorithmWithShift(MatrixNN& H, double tolerance = 1e-20, int max_iterations = 1000) {
    size_t n = H.getSize();
    MatrixNN current_H = H;

    for (int iter = 0; iter < max_iterations; iter++) {
        bool converged = true;
        for (size_t i = 1; i < n; i++) {
            if (std::fabs(current_H[i][i - 1]) > tolerance) {
                converged = false;
                break;
            }
        }

        if (converged) {
            std::cout << "Сходимость достигнута на итерации " << iter + 1 << std::endl;
            break;
        }

        double shift = current_H[n - 1][n - 1];

        MatrixNN shifted_H = current_H;
        for (size_t i = 0; i < n; i++) {
            shifted_H[i][i] -= shift;
        }

        MatrixNN Q = MatrixNN::identity(n);
        MatrixNN R = shifted_H;
        qrDecompositionHessenberg(shifted_H, Q, R);
        // shifted_H.print();

        MatrixNN H_new(n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double sum = 0.0;
                for (size_t k = 0; k < n; k++) {
                    sum += R[i][k] * Q[k][j];
                }
                H_new[i][j] = sum;
            }
            H_new[i][i] += shift;
        }

        current_H = H_new;
    }

    std::vector<double> eigenvalues{};
    for (size_t i = 0; i < n; i++) {
        eigenvalues.push_back(current_H[i][i]);
    }

    std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<double>());

    return eigenvalues;
}

#endif //UNTITLED3_EIGENQR_H
