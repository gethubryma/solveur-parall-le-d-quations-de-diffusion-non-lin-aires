#ifndef METHODES_HH
#define METHODES_HH

#include <vector>

double computeKappa(double u, double k0, double q);
//schéma explicite
void explicitStep(std::vector<double>& u,
                  double dt, double dx,
                  double sigma, double k0, double q,
                  double beta);

//Schéma implicite (linéarisé) 
void implicitStep(std::vector<double>& u,
                  double dt, double dx,
                  double sigma, double k0, double q,
                  double beta);

// Méthode de Newton 
int newtonSolve(std::vector<double>& u,
                double dx, double sigma,
                double k0, double q, double beta,
                double tol, int maxIter);
                
// pour jacobien
std::vector<double> computeResidualNewton(const std::vector<double>& u,
                                          double dx, double sigma,
                                          double k0, double q, double beta);

void computeJacobianNewton(const std::vector<double>& u,
                           double dx, double sigma,
                           double k0, double q, double beta,
                           std::vector<double>& diag,
                           std::vector<double>& upper,
                           std::vector<double>& lower);

// Méthode de Newton - Hypre
int newtonSolveHypre(std::vector<double>& u,
                     double dx, double sigma,
                     double k0, double q, double beta,
                     double tol, int maxIter);

#endif // METHODES_HH
