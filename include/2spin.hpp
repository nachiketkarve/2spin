#ifndef HEADER_2SPIN
#define HEADER_2SPIN

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

//Convert a pair of integers to a single integer and vice versa
int tuple2num(int S1, int S2, int m1, int m2);
void num2tuple(int S1, int S2, int num, int &m1, int &m2);

//Spin operators
double SPlus(int S, int m, int n, double hbar);
double SMinus(int S, int m, int n, double hbar);
std::complex<double> Sx(int S, int m, int n, double hbar);
std::complex<double> Sy(int S, int m, int n, double hbar);
std::complex<double> Sz(int S, int m, int n, double hbar);

//Hamiltonian and other operators
void constructHamiltonian(Matrix &H, int S1, int S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double hbar);
void constructS1zS2z(Matrix &S1zS2z, int S1, int S2, double hbar);
void Update(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt);

//Unitary evolution
void UpdateUnitary(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt);
void UpdateUnitary4thOrder(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt);

//Spectral function
double spectralFunction(const VectorReal &eigVals, const Matrix &eigVecs, const Matrix &V, double omega, double sigma, double beta, double hbar);


#endif
