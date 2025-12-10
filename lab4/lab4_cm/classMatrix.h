#ifndef UNTITLED3_CLASSMATRIX_H
#define UNTITLED3_CLASSMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <utility>

class MatrixNN{
    size_t size;
    std::vector<double> data;

public:

    MatrixNN(): size(0), data() {};

    explicit MatrixNN(const size_t& dim): size(dim){
        data = std::vector<double>(dim*dim, 0.0);
    }; // матрица n размера, заполненная нулями

    explicit MatrixNN(const std::vector<double>& matr)
    {
        size_t s = matr.size();
        size_t root = static_cast<size_t>(std::sqrt(s));
        if (root * root != s) {
            size = 0;
            data.clear();
            throw std::invalid_argument("Матрица не является квадратной");
        } else {
            size = root;
            data = matr;
        }
    }; // матрица из вектора размера nxn

    // Copy constructor (ДОБАВЛЕНО)
    MatrixNN(const MatrixNN& other) : size(other.size), data(other.data) {}

    // Move constructor
    MatrixNN(MatrixNN&& other) noexcept : size(other.size), data(std::move(other.data)) {
        other.size = 0;
    }

    // Copy assignment (ИСПРАВЛЕНО - убрана лишняя реализация)
    MatrixNN& operator=(const MatrixNN& other) = default;

    // Move assignment
    MatrixNN& operator=(MatrixNN&& other) noexcept {
        if (this != &other) {
            size = other.size;
            data = std::move(other.data);
            other.size = 0;
        }
        return *this;
    }

    [[nodiscard]] size_t getSize() const { return size; };
    [[nodiscard]] const std::vector<double>& getData() const { return data; };

    static MatrixNN identity(size_t n) {
        MatrixNN result(n);
        for (size_t i = 0; i < n; i++) {
            result.data[i * n + i] = 1.0;
        }
        return result;
    }

    static MatrixNN matrixFromFile(const std::string& filename)
    {
        MatrixNN result;

        std::ifstream in(filename);
        if (!in.is_open())
        {
            throw std::runtime_error("Не удалось открыть файл: " + filename);
        }

        size_t n;
        in >> n;

        if (!in || n == 0)
        {
            throw std::invalid_argument("Неправильный размер матрицы в файле");
        }

        result = MatrixNN(n);

        for (size_t i = 0; i < n * n; i++)
        {
            in >> result.data[i];
            if (!in)
            {
                throw std::runtime_error("Ошибка чтения элемента " + std::to_string(i));
            }
        }

        return result;
    }

    void print() const {
        for(size_t i = 0; i < size; ++i) {
            for(size_t j = 0; j < size; ++j) {
                std::cout << data[i*size + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    double* operator[](size_t row) {
        if (row >= size || size == 0) {
            throw std::out_of_range("Индекс строки выходит за пределы диапазона");
        }
        return &data[row * size];
    }

    const double* operator[](size_t row) const {
        if (row >= size || size == 0) {
            throw std::out_of_range("Индекс строки выходит за пределы диапазона");
        }
        return &data[row * size];
    }

    /**
     * Решение системы линейных уравнений Ax = b с использованием QR-разложения
     * @param b - вектор правой части системы
     * @return вектор x - решение системы
     */
    std::vector<double> solveQR(const std::vector<double>& b) const {
        if (size == 0) {
            throw std::runtime_error("Матрица пуста");
        }
        if (b.size() != size) {
            throw std::invalid_argument("Размер вектора b не совпадает с размером матрицы");
        }

        // Выполняем QR-разложение
        auto [Q, R] = qrDecomposition();

        // Вычисляем Q^T * b (ИСПРАВЛЕНО)
        std::vector<double> qtb(size, 0.0);
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                qtb[i] += Q[i][j] * b[j]; // Правильное вычисление Q^T * b
            }
        }

        // Решаем систему Rx = Q^T * b обратной подстановкой
        return backSubstitution(R, qtb);
    }

    /**
     * Проверка ортогональности матрицы Q (Q^T * Q = I)
     */
    bool isOrthogonal(const MatrixNN& Q, double tolerance = 1e-10) const {
        size_t n = Q.getSize();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                double dot = 0.0;
                for (size_t k = 0; k < n; ++k) {
                    dot += Q[k][i] * Q[k][j]; // Q^T[i,k] * Q[k,j]
                }
                double expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(dot - expected) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Проверка верхней треугольности матрицы R
     */
    bool isUpperTriangular(const MatrixNN& R, double tolerance = 1e-10) const {
        for (size_t i = 1; i < R.getSize(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (std::abs(R[i][j]) > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

private:
    /**
     * Обратная подстановка для верхней треугольной матрицы
     * @param U - верхняя треугольная матрица
     * @param b - вектор правой части
     * @return решение системы Ux = b
     */
    std::vector<double> backSubstitution(const MatrixNN& U, const std::vector<double>& b) const {
        std::vector<double> x(size, 0.0);

        for (int i = static_cast<int>(size) - 1; i >= 0; --i) {
            // Проверка на вырожденность матрицы (перемещена в начало)
            if (std::abs(U[i][i]) < 1e-12) {
                throw std::runtime_error("Матрица вырождена или близка к вырожденной");
            }

            double sum = 0.0;
            for (size_t j = i + 1; j < size; ++j) {
                sum += U[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / U[i][i];
        }

        return x;
    }

    /**
     * QR-разложение матрицы методом Хаусхолдера
     * @return пара матриц (Q, R)
     */
    std::pair<MatrixNN, MatrixNN> qrDecomposition() const {
        MatrixNN Q = MatrixNN::identity(size);
        MatrixNN R(*this); // Используем copy constructor вместо присваивания

        for (size_t k = 0; k < size - 1; ++k) {
            // Вычисляем вектор x (столбец k, начиная с диагонали)
            std::vector<double> x(size - k);
            for (size_t i = k; i < size; ++i) {
                x[i - k] = R[i][k];
            }

            // Вычисляем норму x
            double norm_x = 0.0;
            for (double val : x) {
                norm_x += val * val;
            }
            norm_x = std::sqrt(norm_x);

            // Пропускаем преобразование если норма слишком мала
            if (norm_x < 1e-12) {
                continue;
            }

            // Создаем вектор v
            std::vector<double> v = x;
            if (x[0] >= 0) {
                v[0] += norm_x;
            } else {
                v[0] -= norm_x;
            }

            // Вычисляем норму v
            double norm_v = 0.0;
            for (double val : v) {
                norm_v += val * val;
            }
            norm_v = std::sqrt(norm_v);

            // Проверка на нулевую норму (добавлена)
            if (norm_v < 1e-12) {
                continue;
            }

            // Нормализуем v
            for (double& val : v) {
                val /= norm_v;
            }

            // Обновляем матрицу R
            for (size_t j = k; j < size; ++j) {
                // Вычисляем скалярное произведение v^T * столбец j
                double dot = 0.0;
                for (size_t i = 0; i < v.size(); ++i) {
                    dot += v[i] * R[k + i][j];
                }

                // Вычитаем 2 * v * (v^T * столбец j)
                for (size_t i = 0; i < v.size(); ++i) {
                    R[k + i][j] -= 2.0 * v[i] * dot;
                }
            }

            // Обновляем матрицу Q
            for (size_t j = 0; j < size; ++j) {
                // Вычисляем скалярное произведение v^T * столбец j Q
                double dot = 0.0;
                for (size_t i = 0; i < v.size(); ++i) {
                    dot += v[i] * Q[k + i][j];
                }

                // Вычитаем 2 * v * (v^T * столбец j Q)
                for (size_t i = 0; i < v.size(); ++i) {
                    Q[k + i][j] -= 2.0 * v[i] * dot;
                }
            }
        }

        return std::make_pair(Q, R); // Используем make_pair явно
    }
};

#endif //UNTITLED3_CLASSMATRIX_H