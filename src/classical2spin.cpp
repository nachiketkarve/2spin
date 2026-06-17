#include "classical2spin.hpp"

double Hamiltonian0(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az)
{
    double Energy = -Jx * S1(0) * S2(0) - Jy * S1(1) * S2(1) - Jz * S1(2) * S2(2) + 0.5 * Ax * (S1(0) * S1(0) + S2(0) * S2(0)) + 0.5 * Ay * (S1(1) * S1(1) + S2(1) * S2(1)) + 0.5 * Az * (S1(2) * S1(2) + S2(2) * S2(2));
    return Energy;
}

double Hamiltonian(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t)
{
    double Energy = Hamiltonian0(S1, S2, Jx, Jy, Jz, Ax, Ay, Az) + lambda * std::sin(mu * t) * std::sin(mu * t) * S1(2) * S2(2);
    return Energy;
}

void EOM(VectorReal &dS1, VectorReal &dS2, VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t)
{
    dS1(0) = (-Jy * S2(1) + Ay * S1(1)) * S1(2) - (-Jz * S2(2) + Az * S1(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S2(2)) * S1(1);
    dS1(1) = (-Jz * S2(2) + Az * S1(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S2(2)) * S1(0) - (-Jx * S2(0) + Ax * S1(0)) * S1(2);
    dS1(2) = (-Jx * S2(0) + Ax * S1(0)) * S1(1) - (-Jy * S2(1) + Ay * S1(1)) * S1(0);

    dS2(0) = (-Jy * S1(1) + Ay * S2(1)) * S2(2) - (-Jz * S1(2) + Az * S2(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S1(2)) * S2(1);
    dS2(1) = (-Jz * S1(2) + Az * S2(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S1(2)) * S2(0) - (-Jx * S1(0) + Ax * S2(0)) * S2(2);
    dS2(2) = (-Jx * S1(0) + Ax * S2(0)) * S2(1) - (-Jy * S1(1) + Ay * S2(1)) * S2(0);
}

VectorReal crossProduct(VectorReal &a, VectorReal &b)
{
    VectorReal c(3);
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

VectorReal rotate(VectorReal &S, VectorReal &B, double angle)
{
    if (B.norm() < TOL)
    {
        return S;
    }
    VectorReal S_rot(3);
    VectorReal Bhat = B / B.norm();
    S_rot = std::cos(angle) * S + std::sin(angle) * crossProduct(Bhat, S) + (1.0 - std::cos(angle)) * (Bhat.dot(S)) * Bhat;
    return S_rot;
}

void Field(VectorReal &B1, VectorReal &B2, VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t)
{
    B1(0) = - Jx * S2(0) + Ax * S1(0);
    B1(1) = - Jy * S2(1) + Ay * S1(1);
    B1(2) = - Jz * S2(2) + Az * S1(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S2(2);

    B2(0) = - Jx * S1(0) + Ax * S2(0);
    B2(1) = - Jy * S1(1) + Ay * S2(1);
    B2(2) = - Jz * S1(2) + Az * S2(2) + lambda * std::sin(mu * t) * std::sin(mu * t) * S1(2);
}

void Update(VectorReal &S1, VectorReal &S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double lambda, double mu, double t, double dt)
{
    VectorReal B1(3), B2(3);
    VectorReal S1mid(3), S2mid(3);
    Field(B1, B2, S1, S2, Jx, Jy, Jz, Ax, Ay, Az, lambda, mu, t);
    S1mid = rotate(S1, B1, B1.norm() * dt / 2.0);
    S2mid = rotate(S2, B2, B2.norm() * dt / 2.0);
    Field(B1, B2, S1mid, S2mid, Jx, Jy, Jz, Ax, Ay, Az, lambda, mu, t + dt / 2.0);
    S1 = rotate(S1, B1, B1.norm() * dt);
    S2 = rotate(S2, B2, B2.norm() * dt);
}