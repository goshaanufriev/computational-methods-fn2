#pragma once

double* simIter(double*& A, double*& b, const int& n)
{
	double* E = idMat(n);
	double k = 0.00001;

	double* A1 = kMat(A, n, n, k);
	double* C = subMat(E, A1, n, n);
	double* y = kMat(b, n, 1, k);
	double* x0 = new double[n];
	for (auto i = 0; i < n; ++i)
	{
		x0[i] = 0;
	}
	double normC = cubeNorm(C, n, n);
	std::cout << normC << "\n";
	double eps = 1e-4;
	double* cx = copyMat(x0, n, 1);
	multMat(C, cx, n, n, n, 1);
	double* x = addMat(cx, y, n, 1);
	double* err = subMat(x, x0, n, 1);
	double errNorm = cubeNorm(err, n, 1);
	int i = 1;
	std::cout << "rho_0 = " << errNorm << "\n";
	while (errNorm > (1 - normC) * eps / normC) {
		x0 = copyMat(x, n, 1);
		cx = copyMat(x0, n, 1);
		multMat(C, cx, n, n, n, 1);
		x = addMat(cx, y, n, 1);
		err = subMat(x, x0, n, 1);
		errNorm = cubeNorm(err, n, 1);
		i++;
	}
	std::cout << i << " " << errNorm << "\n";
	return x;
}

double* jacobi(double*& A, double*& b, const int& n)
{
	double* C = new double[n * n];
	double* y = new double[n];

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

	double* CU = new double[n * n];
	for (auto i = 0; i < n; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			if (j > i)
				CU[i * n + j] = C[i * n + j];
			else
				CU[i * n + j] = 0;
		}
	}

	double eps = 1e-4;
	double normC = cubeNorm(C, n, n);
	double normCU = cubeNorm(CU, n, n);
	std::cout << normC << "\n";

	double* cx = copyMat(x0, n, n);
	multMat(C, cx, n, n, n, 1);
	double* x = addMat(cx, y, n, 1);
	double* err = subMat(x, x0, n, 1);
	double errNorm = cubeNorm(err, n, 1);
	/*double* xk = copyMat(x, n, 1);
	multMat(A, xk, n, n, n, 1);
	double* err = subMat(b, xk, n, 1);
	double errNorm = cubeNorm(err, n);*/
	std::cout << "rho_0 = " << errNorm << "\n";
	int k = 1;
	//while (errNorm > (1 - normC) * eps / normC)
	//while (errNorm > (1 - normC) * eps / normCU)
	while (errNorm > eps)
	//while(k < 13)
	{
		x0 = copyMat(x, n, 1);
		cx = copyMat(x0, n, 1);
		multMat(C, cx, n, n, n, 1);
		x = addMat(cx, y, n, 1);
		err = subMat(x, x0, n, 1);
		errNorm = cubeNorm(err, n, 1);
		/*xk = copyMat(x, n, 1);
		multMat(A, xk, n, n, n, 1);
		err = subMat(b, xk, n, 1);
		errNorm = cubeNorm(err, n);*/
		k += 1;
	}
	std::cout << k << " " << errNorm << "\n";

	return x;
}

double* relax(double*& A, double* b, const int& n, const double& w)
{
	double* L = new double[n * n];
	double* D = new double[n];
	double* U = new double[n * n];

	for (auto i = 0; i < n; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			if (i == j)
			{
				D[i] = A[i * n + j];
				L[i * n + j] = 0;
				U[i * n + j] = 0;
			}
			else if (i > j)
			{
				L[i * n + j] = A[i * n + j];
				U[i * n + j] = 0;
			}
			else
			{
				U[i * n + j] = A[i * n + j];
				L[i * n + j] = 0;
			}
		}
	}
	double* B = kMat(L, n, n, w);
	addToDiagMat(D, B, n);
	double* invB = invMat(n, B);
	double* D1 = kMat(D, n, 1, 1 - w);
	double* C = kMat(U, n, n, -w);
	addToDiagMat(D1, C, n);
	multMat(invB, C, n, n, n, n);
	double* f = copyMat(b, n, 1);
	multMat(invB, f, n, n, n, 1);
	double* y = kMat(f, n, 1, w);

	double* x0 = new double[n];
	for (auto i = 0; i < n; ++i)
	{
		x0[i] = 0;
	}

	double* CU = new double[n * n];
	for (auto i = 0; i < n; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			if (j > i)
				CU[i * n + j] = C[i * n + j];
			else
				CU[i * n + j] = 0;
		}
	}

	double eps = 1e-4;
	double normC = cubeNorm(C, n, n);
	double normCU = cubeNorm(CU, n, n);
	std::cout << "C: " << normC << "\n";
	std::cout << "CU: " << normCU << "\n";

	double* cx = copyMat(x0, n, 1);
	multMat(C, cx, n, n, n, 1);
	double* x = addMat(cx, y, n, 1);
	double* err = subMat(x, x0, n, 1);
	double errNorm = cubeNorm(err, n, 1);
	/*double* xk = copyMat(x, n, 1);
	multMat(A, xk, n, n, n, 1);
	double* err = subMat(b, xk, n, 1);
	double errNorm = cubeNorm(err, n);*/
	std::cout << "rho_0 = " << errNorm << "\n";
	int i = 1;
	//while (errNorm > (1 - normC) * eps / normCU)
	while (errNorm > (1 - normC) * eps / normC)
	//while(errNorm > eps)
	//while(i < 17)
	{
		x0 = copyMat(x, n, 1);
		cx = copyMat(x0, n, 1);
		multMat(C, cx, n, n, n, 1);
		x = addMat(cx, y, n, 1);
		err = subMat(x, x0, n, 1);
		errNorm = cubeNorm(err, n, 1);
		/*xk = copyMat(x, n, 1);
		multMat(A, xk, n, n, n, 1);
		err = subMat(b, xk, n, 1);
		errNorm = cubeNorm(err, n);*/
		++i;
	}
	std::cout << errNorm << " " << i << "\n";

	return x;
}