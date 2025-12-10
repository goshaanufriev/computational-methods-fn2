#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include "classMatrix.h"
#include "eigenqr.h"
#include "inverse_iteration.h"
#include "combinated.h"

// нахождение собственных значений с помощью QR-разложения
void test1()
{
    MatrixNN A = MatrixNN::matrixFromFile("mat1.txt");
    //
    if (A.getSize() == 0) {
        std::cout << "Ошибка: не удалось загрузить матрицу из файла" << std::endl;
        //        return 1;
    }

    std::cout << "\nИсходная матрица:" << std::endl;
    A.print();
    auto Q = MatrixNN::identity(4);

    double trace_A = 0.0;
    for (size_t i = 0; i < A.getSize(); i++) {
        trace_A += A[i][i];
    }

    toHessenberg(A);
    std::cout << "\nМатрица в форме Хессенберга:" << std::endl;
    A.print();

    std::vector<double> eigenvalues = qrAlgorithmWithShift(A);


    std::cout << "\nСобственные значения:" << std::endl;
    for (size_t i = 0; i < eigenvalues.size(); i++) {
        std::cout << "λ" << i + 1 << " = " << std::setprecision(8) << eigenvalues[i] << std::endl;
    }
}

// нахождение собственных значений с помощью QR-разложения + нахождение собственных векторов методом обратных итераций
void test2()
{
    try {
        std::cout << "=== Метод обратной итерации для нахождения собственных векторов ===\n" << std::endl;

        // 1. Считываем матрицу из файла mat.txt
        std::cout << "1. Чтение матрицы из файла mat.txt..." << std::endl;
        MatrixNN A = MatrixNN::matrixFromFile("mat1.txt");

        std::cout << "Размер матрицы: " << A.getSize() << "x" << A.getSize() << std::endl;
        std::cout << "Исходная матрица A:" << std::endl;
        A.print();

        // 2. Приводим матрицу к форме Хессенберга
        std::cout << "\n2. Приведение матрицы к форме Хессенберга..." << std::endl;
        MatrixNN H = A; // Создаем копию матрицы
        toHessenberg(H);
        std::cout << "Матрица в форме Хессенберга:" << std::endl;
        H.print();

        // 3. Находим собственные значения QR-алгоритмом со сдвигом
        std::cout << "\n3. Поиск собственных значений QR-алгоритмом..." << std::endl;
        std::vector<double> eigenvalues = qrAlgorithmWithShift(H);

        std::cout << "Найденные собственные значения: ";
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            std::cout << eigenvalues[i];
            if (i < eigenvalues.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;

        // 4. Находим собственные векторы методом обратной итерации
        std::cout << "\n4. Поиск собственных векторов методом обратной итерации..." << std::endl;
        MatrixNN eigenvectors = inverseIterationMultiple(A, eigenvalues);

        // 5. Выводим результаты
        std::cout << "\n5. Результаты:" << std::endl;
        std::cout << "Собственные векторы (по столбцам):" << std::endl;
        eigenvectors.print();

        // 6. Проверяем точность
        std::cout << "\n6. Проверка точности:" << std::endl;
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            std::vector<double> eigenvector(A.getSize());
            for (size_t j = 0; j < A.getSize(); ++j) {
                eigenvector[j] = eigenvectors[j][i];
            }
            std::cout << "Для λ = " << eigenvalues[i] << ": ";
            verifyEigenPair(A, eigenvector, eigenvalues[i]);
        }

        std::cout << "\n=== Вычисления завершены успешно ===" << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "\n!!! ОШИБКА: " << e.what() << std::endl;
    }
}

//модификацию с использованием отношения Рэлея
void test3()
{
    try {
        std::cout << "=== Метод обратных итераций с отношением Рэлея ===\n" << std::endl;

        MatrixNN A = MatrixNN::matrixFromFile("mat1.txt");

        std::cout << "Исходная матрица A:" << std::endl;
        A.print();
        std::cout << "\n" << std::endl;

        // Находим одну собственную пару методом Рэлея
        auto [eigenvalue, eigenvector] = rayleighQuotientIteration(A);

        // Выводим результаты
        std::cout << "\n" << std::endl;
        std::cout << "Собственное значение: λ = " << std::setprecision(12) << eigenvalue << std::endl;
        std::cout << "Собственный вектор: ";
        for (double val : eigenvector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "\nПроверка точности:" << std::endl;
        verifyEigenPair(A, eigenvalue, eigenvector);


    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }
}

int main() {
    setlocale(LC_ALL, "");

    test2();

    return 0;
}
