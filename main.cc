#include <iostream>
#include <vector>
#include <cmath>
#include "methodes.hh"

void printSolution(const std::vector<double>& u, double dx) {
    for (size_t i = 0; i < u.size(); ++i) {
        std::cout << i * dx << " " << u[i] << "\n";
    }
}

int main() {
    int N = 50;
    double dx = 1.0 / N;
    double sigma = 0.1, k0 = 0.01, q = 1.0, beta = 1.0;
    double gamma = 0.1;
    double dt = gamma * 2.0 / (4.0 * sigma * std::pow(2.0, 3)
                              + 4.0 * k0 * std::pow(2.0, q) / (dx * dx));
    int nsteps = 1000;

    {
        std::vector<double> u_explicit(N+1, 1.0);
        for (int n = 0; n < nsteps; ++n) {
            explicitStep(u_explicit, dt, dx, sigma, k0, q, beta);
        }
        std::cout << "Solution finale u(x) (Schéma explicite) \n";
        // printSolution(u_explicit, dx);
        //printAmplificationMatrixExplicit(u_explicit, dt, dx, sigma, k0, q);
    }

    {
        std::vector<double> u_implicit(N+1, 1.0);
        for (int n = 0; n < nsteps; ++n) {
            implicitStep(u_implicit, dt, dx, sigma, k0, q, beta);
        }
        std::cout << "\nSolution finale u(x) (Schéma implicite linéarisé) \n";
        //printSolution(u_implicit, dx);
        //printImplicitMatrix(u_implicit, dt, dx, sigma, k0, q, beta);
    }

    {
        std::vector<double> u_newton(N+1, 1.0);
        double tol = 1.0e-10;
        int maxIter = 50;
        
        std::cout << "\nSolution finale u(x) (Méthode de Newton) :\n";
        newtonSolve(u_newton, dx, sigma, k0, q, beta, tol, maxIter);
        //printSolution(u_newton, dx);
    }

    return 0;
}
