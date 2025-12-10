#ifndef UNTITLED3_COMBINATED_H
#define UNTITLED3_COMBINATED_H

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>
#include <iomanip>
#include <utility>

// Вычисление отношения Рэлея для вектора
double rayleighQuotient(const MatrixNN& A, const std::vector<double>& x) {
    size_t n = A.getSize();
    if (x.size() != n) {
        throw std::invalid_argument("Размер вектора не совпадает с размером матрицы");
    }

    // Вычисляем Ax
    std::vector<double> Ax(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    // Вычисляем (Ax, x) и (x, x)
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < n; ++i) {
        numerator += Ax[i] * x[i];
        denominator += x[i] * x[i];
    }

    if (std::abs(denominator) < 1e-14) {
        throw std::runtime_error("Нулевой вектор в отношении Рэлея");
    }

    return numerator / denominator;
}

//Модифицированный метод обратных итераций с использованием отношения Рэлея
std::pair<double, std::vector<double>> rayleighQuotientIteration(
        const MatrixNN& A,
        size_t max_iterations = 1000,
        double tolerance = 1e-12) {

    size_t n = A.getSize();
    if (n == 0) {
        throw std::invalid_argument("Матрица пуста");
    }

    // Генератор случайных чисел для начального приближения
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    // Начальное приближение - случайный нормированный вектор
    std::vector<double> x(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = dis(gen);
    }

    // Нормализация начального вектора
    double norm = 0.0;
    for (double val : x) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    for (double& val : x) {
        val /= norm;
    }

    // Начальное приближение для собственного значения через отношение Рэлея
    double lambda = rayleighQuotient(A, x);

    std::cout << "Начальное приближение: λ = " << std::setprecision(12) << lambda << std::endl;

    // Итерационный процесс
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        // Создаем матрицу (A - λI)
        MatrixNN A_minus_lambdaI(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                A_minus_lambdaI[i][j] = A[i][j];
            }
            A_minus_lambdaI[i][i] -= lambda;
        }

        // Решаем систему (A - λI)y = x
        std::vector<double> y;
        try {
            y = A_minus_lambdaI.solveQR(x);
        } catch (const std::exception& e) {
            // Если матрица вырождена, добавляем небольшое возмущение
            MatrixNN perturbed = A_minus_lambdaI;
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

        if (new_norm < 1e-14) {
            // Если норма слишком мала, перезапускаем со случайным вектором
            for (size_t i = 0; i < n; ++i) {
                x[i] = dis(gen);
            }
            norm = 0.0;
            for (double val : x) {
                norm += val * val;
            }
            norm = std::sqrt(norm);
            for (double& val : x) {
                val /= norm;
            }
            lambda = rayleighQuotient(A, x);
            continue;
        }

        // Сохраняем старый вектор для проверки сходимости
        std::vector<double> x_old = x;

        // Обновляем вектор
        for (size_t i = 0; i < n; ++i) {
            x[i] = y[i] / new_norm;
        }

        // Обновляем собственное значение через отношение Рэлея
        double new_lambda = rayleighQuotient(A, x);

        // Проверяем сходимость по изменению собственного значения и вектора
        double delta_lambda = std::abs(new_lambda - lambda);

        // Проверяем изменение вектора
        double delta_x = 0.0;
        for (size_t i = 0; i < n; ++i) {
            delta_x += std::abs(x[i] - x_old[i]);
        }

        lambda = new_lambda;

        if (delta_lambda < tolerance && delta_x < tolerance) {
            std::cout << "Метод Рэлея сошелся за " << iter + 1 << " итераций" << std::endl;
            std::cout << "Найденное собственное значение: " << std::setprecision(12) << lambda << std::endl;
            break;
        }

        if (iter == max_iterations - 1) {
            std::cout << "Предупреждение: достигнуто максимальное количество итераций" << std::endl;
            std::cout << "Текущее приближение: λ = " << std::setprecision(12) << lambda << std::endl;
        }
    }

    return {lambda, x};
}

// Проверка точности собственной пары
void verifyEigenPair(const MatrixNN& A, double eigenvalue, const std::vector<double>& eigenvector,
                     double tolerance = 1e-10) {
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

    // Вычисляем относительную невязку
    double residual = 0.0;
    double norm_Ax = 0.0;
    for (size_t i = 0; i < n; ++i) {
        residual += (Ax[i] - lambda_x[i]) * (Ax[i] - lambda_x[i]);
        norm_Ax += Ax[i] * Ax[i];
    }
    residual = std::sqrt(residual);
    norm_Ax = std::sqrt(norm_Ax);

    double relative_residual = residual / (norm_Ax + std::abs(eigenvalue));

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Относительная невязка: " << relative_residual;
    if (relative_residual < tolerance) {
        std::cout << " (OK)" << std::endl;
    } else {
        std::cout << " (ВНИМАНИЕ: большая невязка)" << std::endl;
    }
    std::cout << std::defaultfloat;
}

#endif //UNTITLED3_COMBINATED_H