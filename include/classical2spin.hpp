#ifndef HEADER_CLASSICAL2SPIN
#define HEADER_CLASSICAL2SPIN

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include "eigenClasses.hpp"

const double TOL = 1e-12;


// Classical 2-spin system
double Hamiltonian0(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az);
double Hamiltonian(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t);

// Equations of motion
void EOM(VectorReal &dS1, VectorReal &dS2, VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t);

// Vector operations
VectorReal crossProduct(VectorReal &a, VectorReal &b);
VectorReal rotate(VectorReal &S, VectorReal &B, double angle);

// Magnetic field
void Field(VectorReal &B1, VectorReal &B2, VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t);

// Update the state of the system
void Update(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t, double dt);












#endif