#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

double PI = 2 * acos(0.0);

std::vector<double> uniformGrid(const double& a, const double& b, const int& n)
{
    std::vector<double> grid(n+1);
    double h = (b - a) / n;
    for (auto i = 0; i <= n; ++i)
    {
        grid[i] = a + i * h;
    }
    return grid;
}

std::vector<double> chebGrid(const double& a, const double& b, const int& n)
{
    std::vector<double> grid(n + 1);
    for (auto i = 0; i <= n; ++i)
    {
        grid[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * PI / (2 * (n + 1)));
    }
    return grid;
}

std::vector<double> yValues(const std::vector<double>& grid)
{
    std::vector<double> y(grid.size());
    for (auto i = 0; i < grid.size(); ++i)
    {
        //y[i] = grid[i] * grid[i]; // x^2
        //y[i] = 1 / (1 + grid[i] * grid[i]);
        y[i] = exp(grid[i]);
    }
    return y;
}

double Lagranj(const double& x, const std::vector<double>& grid, const std::vector<double>& values)
{
    double L = 0;

    for (auto i = 0; i < grid.size(); ++i)
    {
        double c = 1;
        for (auto j = 0; j < grid.size(); ++j)
        {
            if (i != j)
                c *= (x - grid[j]) / (grid[i] - grid[j]);
        }
        L += c * values[i];
    }

    return L;
}

std::vector<double> tridiag(const int& n, const std::vector<double>& h, const std::vector<double>& f)
{
    std::vector<double> t(n - 1);
    std::vector<double> alpha(n);
    std::vector<double> beta(n);
    alpha[1] = -h[1] / (h[0] + h[1]);
    beta[1] = f[1] / (h[0] + h[1]);
    beta[n - 1] = (f[n - 1] - h[n - 2] * beta[n - 2]) / (h[n - 2] * alpha[n - 2] + 2 * h[n - 2] + 2 * h[n - 1]);
    t[n - 2] = beta[n - 1];
    for (auto i = 0; i < n - 1; ++i)
    {
        double denominator = h[i] * alpha[i] + 2 * h[i] + 2 * h[i + 1];
        alpha[i + 1] = -h[i + 1] / denominator;
        beta[i + 1] = f[i+1]-h[i]*beta[i] / denominator;
    }
    for (auto i = n - 3; i > -1; --i)
    {
        t[i] = alpha[i + 1] * t[i + 1] + beta[i + 1];
    }
    return t;
}

void KV3()
{
    double a = 0;
    double b = 2;
    double h = 0.2;
    double n = (b - a) / h;

    std::vector<double> xgrid = uniformGrid(a, b, n);
    std::vector<double> ygrid = yValues(xgrid);
    double x = 2.2;
    double y = Lagranj(x, xgrid, ygrid);

    double err = std::fabs(y - exp(x));
    std::cout << err << "\n";
}

int main()
{
    setlocale(LC_ALL, "Russian");

    std::ifstream fin("input.txt");
    int n, a, b;
    fin >> n >> a >> b;

    std::vector<double> grid = chebGrid(a, b, n);
    std::vector<double> values = yValues(grid);

    double x = 0.2;
    double y = Lagranj(x, grid, values);

    KV3();

    std::ofstream fout("output.txt");
    fout << y;

    return 0;
}