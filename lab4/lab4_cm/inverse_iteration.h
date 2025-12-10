#ifndef UNTITLED3_INVERSE_ITERATION_H
#define UNTITLED3_INVERSE_ITERATION_H


#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

std::vector<double> inverseIteration(const MatrixNN& A, double eigenvalue_approx,
                                     size_t max_iterations = 1000, double tolerance = 1e-10) {
    size_t n = A.getSize();
    if (n == 0) {
        throw std::invalid_argument("Empty matrix");
    }

    // Создаем матрицу (A - μI), где μ - приближение собственного значения
    MatrixNN A_minus_muI(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A_minus_muI[i][j] = A[i][j];
        }
        A_minus_muI[i][i] -= eigenvalue_approx; // вычитаем μ на диагонали
    }

    // Генератор случайных чисел для начального приближения
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    // Начальное приближение случайным вектором
    std::vector<double> x(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = dis(gen);
    }

    // Нормализуем начальный вектор
    double norm = 0.0;
    for (double val : x) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    for (double& val : x) {
        val /= norm;
    }

    // Итерационный процесс
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        // Решаем систему (A - μI)y = x
        std::vector<double> y;
        try {
            y = A_minus_muI.solveQR(x);
        } catch (const std::exception& e) {
            // Если матрица вырождена, добавляем небольшое возмущение
            MatrixNN perturbed = A_minus_muI;
            for (size_t i = 0; i < n; ++i) {
                perturbed[i][i] += 1e-12;
            }
            y = perturbed.solveQR(x);
        }

        // Нормализуем новый вектор
        double new_norm = 0.0;
        for (double val : y) {
            new_norm += val * val;
        }
        new_norm = std::sqrt(new_norm);

        // Проверяем сходимость
        double diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double normalized_val = y[i] / new_norm;
            diff += std::abs(normalized_val - x[i]);
            x[i] = normalized_val;
        }

        if (diff < tolerance) {
            std::cout << "Inverse iteration: " << iter + 1 << "iterations" << std::endl;
            break;
        }

        if (iter == max_iterations - 1) {
            std::cout << "Предупреждение: достигнуто максимальное количество итераций" << std::endl;
        }
    }

    return x;
}


MatrixNN inverseIterationMultiple(const MatrixNN& A, const std::vector<double>& eigenvalues,
                                  size_t max_iterations = 1000, double tolerance = 1e-10) {
    size_t n = A.getSize();
    if (eigenvalues.size() != n) {
        throw std::invalid_argument("Количество приближений собственных значений должно совпадать с размером матрицы");
    }

    MatrixNN eigenvectors(n);

    for (size_t i = 0; i < n; ++i) {
        std::cout << "Approx eigenvalue" << eigenvalues[i] << "..." << std::endl;

        std::vector<double> eigenvector = inverseIteration(A, eigenvalues[i], max_iterations, tolerance);

        // Сохраняем собственный вектор в i-й столбец матрицы
        for (size_t j = 0; j < n; ++j) {
            eigenvectors[j][i] = eigenvector[j];
        }
    }

    return eigenvectors;
}


void verifyEigenPair(const MatrixNN& A, const std::vector<double>& eigenvector, double eigenvalue,
                     double tolerance = 1e-8) {
    size_t n = A.getSize();
    std::vector<double> Ax(n, 0.0);
    std::vector<double> lambda_x(n, 0.0);

    // Вычисляем Ax
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Ax[i] += A[i][j] * eigenvector[j];
        }
    }

    // Вычисляем λx
    for (size_t i = 0; i < n; ++i) {
        lambda_x[i] = eigenvalue * eigenvector[i];
    }

    // Вычисляем невязку ||Ax - λx||
    double residual = 0.0;
    for (size_t i = 0; i < n; ++i) {
        residual += std::abs(Ax[i] - lambda_x[i]);
    }

    std::cout << "Невязка ||Ax - λx||: " << residual;
    if (residual < tolerance) {
        std::cout << " (OK)" << std::endl;
    } else {
        std::cout << " (ВНИМАНИЕ: большая невязка)" << std::endl;
    }
}

#endif //UNTITLED3_INVERSE_ITERATION_H
