#include "methodes.hh"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"


double computeKappa(double u, double k0, double q) {
    return k0 * std::pow(u, q);
}

// Schéma explicite
void explicitStep(std::vector<double>& u, double dt, double dx,
                  double sigma, double k0, double q, double beta)
{
    int N = static_cast<int>(u.size()) - 1;
    std::vector<double> unew(u.size(), 0.0);

    double u_mirror = u[0];
    double Q0 = (0.0 <= 0.2) ? beta : 0.0;

    double kappa_ip = 0.5 * (computeKappa(u[0], k0, q) + computeKappa(u[1], k0, q));
    double kappa_im = 0.5 * (computeKappa(u[0], k0, q) + computeKappa(u_mirror, k0, q));
    double conduction = (kappa_ip*(u[1] - u[0]) - kappa_im*(u[0] - u_mirror)) / (dx*dx);
    double radiation  = sigma*(std::pow(u[0], 4.0) - 1.0);
    unew[0] = u[0] + dt*(conduction - radiation + Q0);

    for (int i = 1; i < N; ++i) {
        double x  = i * dx;
        double Qi = (x <= 0.2) ? beta : 0.0;

        double kappa_ip = 0.5 * (computeKappa(u[i],   k0, q) + computeKappa(u[i+1], k0, q));
        double kappa_im = 0.5 * (computeKappa(u[i],   k0, q) + computeKappa(u[i-1], k0, q));
        double conduction = (kappa_ip*(u[i+1] - u[i]) - kappa_im*(u[i] - u[i-1])) / (dx*dx);
        double radiation  = sigma*(std::pow(u[i], 4.0) - 1.0);

        unew[i] = u[i] + dt*(conduction - radiation + Qi);
    }

    // condition dirichlet u(N)=1
    unew[N] = 1.0;
    u = unew;
}

//Schéma implicite linéarisé 
static void solveTridiag(std::vector<double>& a,
                         const std::vector<double>& b,
                         std::vector<double>& c,
                         std::vector<double>& d)
{
    //Thomas
    int N = (int)a.size();
    std::vector<double> c_star(N, 0.0), d_star(N, 0.0);

    c_star[0] = b[0] / a[0];
    d_star[0] = d[0] / a[0];

    for(int i = 1; i < N; i++){
        double m = a[i] - c[i]*c_star[i-1];
        if(i < N) c_star[i] = b[i] / m;
        d_star[i] = (d[i] - c[i]*d_star[i-1]) / m;
    }
    d[N-1] = d_star[N-1];
    for(int i = N-2; i >= 0; --i){
        d[i] = d_star[i] - c_star[i]*d[i+1];
    }
}

void implicitStep(std::vector<double>& u, double dt, double dx,
                  double sigma, double k0, double q, double beta)
{
    int N = static_cast<int>(u.size()) - 1;
    std::vector<double> a(N+1, 0.0), b_diag(N+1, 0.0), c(N+1, 0.0), d(N+1, 0.0);

    double kappa_half = 0.5*(computeKappa(u[0], k0, q) + computeKappa(u[1], k0, q));
    a[0]     = 1 + dt*(2*kappa_half/(dx*dx) + 4*sigma*std::pow(u[0],3));
    b_diag[0]= -dt*(2*kappa_half/(dx*dx));
    d[0]     = u[0] + dt*(beta + sigma*(std::pow(u[0],4)*3 + 1));

    for(int i = 1; i < N; ++i) {
        double x  = i*dx;
        double Qi = (x <= 0.2) ? beta : 0.0;

        double kappa_ip = 0.5*(computeKappa(u[i],   k0, q) + computeKappa(u[i+1], k0, q));
        double kappa_im = 0.5*(computeKappa(u[i],   k0, q) + computeKappa(u[i-1], k0, q));

        a[i]     = 1 + dt*((kappa_ip+kappa_im)/(dx*dx) + 4*sigma*std::pow(u[i],3));
        b_diag[i]= -dt*(kappa_ip/(dx*dx));
        c[i]     = -dt*(kappa_im/(dx*dx));
        d[i]     = u[i] + dt*(Qi + sigma*(3*std::pow(u[i],4) + 1));
    }

    a[N] = 1.0;
    d[N] = 1.0;
    solveTridiag(a, b_diag, c, d);
    u = d;
}
// Méthode de Newton 
std::vector<double> computeResidualNewton(const std::vector<double>& u,
                                          double dx, double sigma,
                                          double k0, double q, double beta)
{
    int N = (int)u.size()-1;
    std::vector<double> F(N+1, 0.0);

    double u_mirror = u[1];
    double kappa_im = 0.5*(computeKappa(u[0], k0, q)+computeKappa(u_mirror, k0, q));
    double kappa_ip = 0.5*(computeKappa(u[0], k0, q)+computeKappa(u[1], k0, q));

    double conduction_0 = (kappa_ip*(u[1]-u[0]) - kappa_im*(u[0]-u_mirror))/(dx*dx);
    double radiation_0  = sigma*(std::pow(u[0],4)-1.0);
    double Q0           = (0.0<=0.2)? beta : 0.0;
    F[0] = - conduction_0 + radiation_0 - Q0;

    for(int i=1; i<N; ++i) {
        double x  = i*dx;
        double Qi = (x <= 0.2)? beta : 0.0;

        double kappa_left  = 0.5*(computeKappa(u[i], k0, q)+computeKappa(u[i-1], k0, q));
        double kappa_right = 0.5*(computeKappa(u[i], k0, q)+computeKappa(u[i+1], k0, q));

        double conduction  = (kappa_right*(u[i+1]-u[i]) - kappa_left*(u[i]-u[i-1]))/(dx*dx);
        double radiation   = sigma*(std::pow(u[i],4)-1.0);
        F[i] = - conduction + radiation - Qi;
    }
    F[N] = u[N] - 1.0;
    return F;
}

void computeJacobianNewton(const std::vector<double>& u,
                           double dx, double sigma,
                           double k0, double q, double beta,
                           std::vector<double>& diag,
                           std::vector<double>& upper,
                           std::vector<double>& lower )
{
    int N = (int)u.size()-1;
    diag.assign(N+1, 0.0);
    upper.assign(N+1, 0.0);
    lower.assign(N+1, 0.0);
    {
        double u_mirror = u[1];
        double kappa_im = 0.5*(computeKappa(u[0],k0,q)+computeKappa(u_mirror,k0,q));
        double kappa_ip = 0.5*(computeKappa(u[0],k0,q)+computeKappa(u[1],k0,q));
        double d_conduc = (kappa_im + kappa_ip)/(dx*dx);
        double d_radiat= 4.0*sigma*std::pow(u[0],3);

        diag[0]  = d_conduc + d_radiat;     
        upper[0] = -kappa_ip/(dx*dx);       
        lower[0] = -kappa_im/(dx*dx);       
    }

    for(int i=1; i<N; ++i) {
        double kappa_im = 0.5*(computeKappa(u[i],k0,q)+computeKappa(u[i-1],k0,q));
        double kappa_ip = 0.5*(computeKappa(u[i],k0,q)+computeKappa(u[i+1],k0,q));
        double d_conduc = (kappa_im + kappa_ip)/(dx*dx);
        double d_radiat= 4.0*sigma*std::pow(u[i],3);

        diag[i]  = d_conduc + d_radiat;
        upper[i] = -kappa_ip/(dx*dx);
        lower[i] = -kappa_im/(dx*dx);
    }
    diag[N]  = 1.0;
    upper[N] = 0.0;
    lower[N] = 0.0;
}

int newtonSolve(std::vector<double>& u,
                double dx, double sigma, double k0, double q, double beta,
                double tol, int maxIter)
{
    int N = (int)u.size()-1;
    for(int iter = 0; iter < maxIter; ++iter) {
        std::vector<double> F = computeResidualNewton(u, dx, sigma, k0, q, beta);

        std::vector<double> diag, upper, lower;
        computeJacobianNewton(u, dx, sigma, k0, q, beta, diag, upper, lower);

        for(int i=0; i<=N; ++i) {
            F[i] = -F[i];
        }
        std::vector<double> a = diag;
        std::vector<double> b = upper;
        std::vector<double> c = lower;
        std::vector<double> d = F;

        solveTridiag(a, b, c, d);

        double maxDelta = 0.0;
        for(int i=0; i<=N; ++i) {
            double delta = d[i];
            u[i] += delta;
            double ad = std::fabs(delta);
            if(ad > maxDelta) maxDelta = ad;
        }
        if(maxDelta < tol) {
            std::cout << "Newton converge en " << iter+1
                      << " itérations, ||delta||inf=" << maxDelta << "\n";
            return iter+1;
        }
    }
    std::cout << "Newton n'a pas convergé en " << maxIter << " itérations.\n";
    return maxIter;
}

//Méthode de Newton - Hypre
int newtonSolveHypre(std::vector<double>& u,
                     double dx, double sigma,
                     double k0, double q, double beta,
                     double tol, int maxIter)
{
    int N = (int)u.size()-1;

    for(int iter=0; iter<maxIter; ++iter)
    {
        //Calculer du résidu
        std::vector<double> F = computeResidualNewton(u, dx, sigma, k0, q, beta);

        std::vector<double> diag, upper, lower;
        computeJacobianNewton(u, dx, sigma, k0, q, beta, diag, upper, lower);

        for(int i=0; i<=N; ++i) {
            F[i] = -F[i];
        }

        int ilower = 0;
        int iupper = N; 

        HYPRE_IJMatrix A;
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
        HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A);

        for(int i = ilower; i <= iupper; i++) {
            std::vector<int>    cols;
            std::vector<double> vals;

            cols.push_back(i);
            vals.push_back(diag[i]);

            if(i>0) {
                cols.push_back(i-1);
                vals.push_back(lower[i]);
            }

            if(i<N) {
                cols.push_back(i+1);
                vals.push_back(upper[i]);
            }

            int ncols = (int)cols.size();
            HYPRE_IJMatrixSetValues(A, 1, &ncols, &i, cols.data(), vals.data());
        }
        HYPRE_IJMatrixAssemble(A);

        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);

        HYPRE_IJVector b_;
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b_);
        HYPRE_IJVectorSetObjectType(b_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_);

        for(int i = ilower; i <= iupper; i++) {
            double val = F[i];
            HYPRE_IJVectorSetValues(b_, 1, &i, &val);
        }
        HYPRE_IJVectorAssemble(b_);
        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(b_, (void**)&par_b);

        HYPRE_IJVector x_;
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x_);
        HYPRE_IJVectorSetObjectType(x_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x_);

        for(int i = ilower; i <= iupper; i++){
            double initVal = 0.0;
            HYPRE_IJVectorSetValues(x_, 1, &i, &initVal);
        }
        HYPRE_IJVectorAssemble(x_);
        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(x_, (void**)&par_x);

        HYPRE_Solver solver;
        HYPRE_BoomerAMGCreate(&solver);

        HYPRE_BoomerAMGSetCoarsenType(solver, 8);    
        HYPRE_BoomerAMGSetInterpType(solver, 6);     
        HYPRE_BoomerAMGSetPrintLevel(solver, 1);     
        HYPRE_BoomerAMGSetRelaxType(solver, 6);    

        HYPRE_BoomerAMGSetTol(solver, 1e-14);
        HYPRE_BoomerAMGSetMaxIter(solver, 50);

        HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

        std::vector<double> delta(N+1, 0.0);
        {
            std::vector<int> indices(N+1);
            for(int i=0; i<=N; ++i) {
                indices[i] = i;
            }
            HYPRE_IJVectorGetValues(x_, N+1, indices.data(), delta.data());
        }

        //nettoyer
        HYPRE_BoomerAMGDestroy(solver);
        HYPRE_IJMatrixDestroy(A);
        HYPRE_IJVectorDestroy(b_);
        HYPRE_IJVectorDestroy(x_);

        //Mise à jour de u
        double maxDelta = 0.0;
        for(int i=0; i<=N; ++i) {
            u[i] += delta[i];
            double ad = std::fabs(delta[i]);
            if(ad > maxDelta) maxDelta = ad;
        }

        if(maxDelta < tol) {
            std::cout << "Newton converge en " << iter+1
                      << " itérations, ||delta||inf=" << maxDelta << "\n";
            return iter+1;
        }
    }

    std::cout<<"Newton n'a pas convergé en "<<maxIter<<" itérations.\n";
    return maxIter;
}
