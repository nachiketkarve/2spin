#include "2spin.hpp"

int tuple2num(int S1, int S2, int m1, int m2)
{
    // Convert the tuple (m1, m2) to a unique integer index
    return (S1 + m1) * (2 * S2 + 1) + (S2 + m2);
}

void num2tuple(int S1, int S2, int num, int &m1, int &m2)
{
    // Convert the integer index back to the tuple (m1, m2)
    m1 = num / (2 * S2 + 1) - S1;
    m2 = num % (2 * S2 + 1) - S2;
}

double SPlus(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S+ |n>
    if (m != n + 1)
    {
        return 0.0;
    }
    return hbar * std::sqrt(double(S * (S + 1) - n * (n + 1)));
}

double SMinus(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S- |n>
    if (m != n - 1)
    {
        return 0.0;
    }
    return hbar * std::sqrt(double(S * (S + 1) - n * (n - 1)));
}

double SPlus2(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S+^2 |n>
    if (m != n + 2)
    {
        return 0.0;
    }
    return hbar * hbar * std::sqrt(double((S * (S + 1) - n * (n + 1)) * (S * (S + 1) - (n + 1) * (n + 2))));
}

double SMinus2(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S-^2 |n>
    if (m != n - 2)
    {
        return 0.0;
    }
    return hbar * hbar * std::sqrt(double((S * (S + 1) - n * (n - 1)) * (S * (S + 1) - (n - 1) * (n - 2))));
}

double SPlusSMinus(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S+ S- |n>
    if (m != n)
    {
        return 0.0;
    }
    return hbar * hbar * (S * (S + 1) - n * (n - 1));
}

double SMinusSPlus(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| S- S+ |n>
    if (m != n)
    {
        return 0.0;
    }
    return hbar * hbar * (S * (S + 1) - n * (n + 1));
}

int delta(int m, int n)
{
    // Kronecker delta function
    return (m == n) ? 1 : 0;
}

std::complex<double> Sx(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| Sx |n>
    return 0.5 * (SPlus(S, m, n, hbar) + SMinus(S, m, n, hbar));
}

std::complex<double> Sy(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| Sy |n>
    return 0.5 * (SPlus(S, m, n, hbar) - SMinus(S, m, n, hbar)) * std::complex<double>(0, -1);
}

std::complex<double> Sz(int S, int m, int n, double hbar)
{
    // Calculate the matrix element <m| Sz |n>
    (void)S;
    if (m != n)
    {
        return 0.0;
    }
    return hbar * m;
}

void constructHamiltonian(Matrix &H, int S1, int S2, double Jx, double Jy, double Jz, double Ax, double Ay, double Az, double hbar)
{
    int dim = (2 * S1 + 1) * (2 * S2 + 1);
    if (H.rows() != dim || H.cols() != dim)
    {
        throw std::invalid_argument("Matrix H must be of size (2*S1+1)*(2*S2+1) x (2*S1+1)*(2*S2+1)");
    }

    H.setZero();
    
    double AP = 0.125 * (Ax + Ay);
    double AM = 0.125 * (Ax - Ay);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < dim; i++)
    {
        int m1, m2;
        num2tuple(S1, S2, i, m1, m2);

        for (int j = 0; j < dim; j++)
        {
            int n1, n2;
            num2tuple(S1, S2, j, n1, n2);
            std::complex<double> S1x = Sx(S1, m1, n1, hbar);
            std::complex<double> S1y = Sy(S1, m1, n1, hbar);
            std::complex<double> S1z = Sz(S1, m1, n1, hbar);
            std::complex<double> S2x = Sx(S2, m2, n2, hbar);
            std::complex<double> S2y = Sy(S2, m2, n2, hbar);
            std::complex<double> S2z = Sz(S2, m2, n2, hbar);
            std::complex<double> S1PS1P = SPlus2(S1, m1, n1, hbar) * delta(m2, n2);
            std::complex<double> S1MS1M = SMinus2(S1, m1, n1, hbar) * delta(m2, n2);
            std::complex<double> S2PS2P = SPlus2(S2, m2, n2, hbar) * delta(m1, n1);
            std::complex<double> S2MS2M = SMinus2(S2, m2, n2, hbar) * delta(m1, n1);
            std::complex<double> S1PS1M = SPlusSMinus(S1, m1, n1, hbar) * delta(m2, n2);
            std::complex<double> S1MS1P = SMinusSPlus(S1, m1, n1, hbar) * delta(m2, n2);
            std::complex<double> S2PS2M = SPlusSMinus(S2, m2, n2, hbar) * delta(m1, n1);
            std::complex<double> S2MS2P = SMinusSPlus(S2, m2, n2, hbar) * delta(m1, n1);
            std::complex<double> S1z2 = S1z * S1z * double(delta(m2, n2));
            std::complex<double> S2z2 = S2z * S2z * double(delta(m1, n1));

            H(i, j) = - Jx * S1x * S2x - Jy * S1y * S2y - Jz * S1z * S2z + AM * (S1PS1P + S2PS2P + S1MS1M + S2MS2M) + AP * (S1PS1M + S1MS1P + S2PS2M + S2MS2P) + 0.5 * Az * (S1z2 + S2z2);
        }
    }

    Matrix Hd = H.adjoint();
    H = 0.5 * (H + Hd); // Ensure the Hamiltonian is Hermitian
}

void constructS1zS2z(Matrix &S1zS2z, int S1, int S2, double hbar)
{
    int dim = (2 * S1 + 1) * (2 * S2 + 1);
    if (S1zS2z.rows() != dim || S1zS2z.cols() != dim)
    {
        throw std::invalid_argument("Matrix S1zS2z must be of size (2*S1+1)*(2*S2+1) x (2*S1+1)*(2*S2+1)");
    }

    S1zS2z.setZero();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < dim; i++)
    {
        int m1, m2;
        num2tuple(S1, S2, i, m1, m2);

        for (int j = 0; j < dim; j++)
        {
            int n1, n2;
            num2tuple(S1, S2, j, n1, n2);
            std::complex<double> S1z = Sz(S1, m1, n1, hbar);
            std::complex<double> S2z = Sz(S2, m2, n2, hbar);

            S1zS2z(i, j) = S1z * S2z;
        }
    }

    Matrix S1zS2z_d = S1zS2z.adjoint();
    S1zS2z = 0.5 * (S1zS2z + S1zS2z_d); // Ensure the matrix is Hermitian
}

void Update(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt)
{
    Matrix K1 = -std::complex<double>(0, 1) * (H + lambda * std::sin(mu * t) * std::sin(mu * t) * S1zS2z) * U / hbar;
    Matrix K2 = -std::complex<double>(0, 1) * (H + lambda * std::sin(mu * (t + 0.5 * dt)) * std::sin(mu * (t + 0.5 * dt)) * S1zS2z) * (U + 0.5 * dt * K1) / hbar;
    Matrix K3 = -std::complex<double>(0, 1) * (H + lambda * std::sin(mu * (t + 0.5 * dt)) * std::sin(mu * (t + 0.5 * dt)) * S1zS2z) * (U + 0.5 * dt * K2) / hbar;
    Matrix K4 = -std::complex<double>(0, 1) * (H + lambda * std::sin(mu * (t + dt)) * std::sin(mu * (t + dt)) * S1zS2z) * (U + dt * K3) / hbar;
    U += (dt / 6.0) * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
}

void UpdateUnitary(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt)
{
    double tm = t + 0.5 * dt;
    Matrix Hmid = H + lambda * std::pow(std::sin(mu * tm), 2) * S1zS2z;

    Eigen::SelfAdjointEigenSolver<Matrix> eig(Hmid);
    VectorReal evals = eig.eigenvalues();
    Matrix evecs = eig.eigenvectors();

    Matrix phase = Matrix::Zero(H.rows(), H.cols());
    for (int i = 0; i < evals.size(); ++i)
    {
        phase(i, i) = std::exp(-std::complex<double>(0, 1) * evals(i) * dt / hbar);
    }

    Matrix Ustep = evecs * phase * evecs.adjoint();
    U = Ustep * U;
}

void UpdateUnitary4thOrder(const Matrix &H, const Matrix &S1zS2z, Matrix &U, double hbar, double lambda, double mu, double t, double dt)
{
    double a1 = 0.25 + std::sqrt(3.0) / 6.0;
    double a2 = 0.25 - std::sqrt(3.0) / 6.0;
    double c1 = 0.5 - std::sqrt(3.0) / 6.0;
    double c2 = 0.5 + std::sqrt(3.0) / 6.0;

    double t1 = t + c1 * dt;
    double t2 = t + c2 * dt;

    Matrix A1 = H + lambda * std::pow(std::sin(mu * t1), 2) * S1zS2z;
    Matrix A2 = H + lambda * std::pow(std::sin(mu * t2), 2) * S1zS2z;

    Matrix K1 = dt * (a2 * A1 + a1 * A2) / hbar;
    Matrix K2 = dt * (a1 * A1 + a2 * A2) / hbar;

    Eigen::SelfAdjointEigenSolver<Matrix> eig1(K1);
    VectorReal evals1 = eig1.eigenvalues();
    Matrix evecs1 = eig1.eigenvectors();

    Matrix phase1 = Matrix::Zero(H.rows(), H.cols());
    for (int i = 0; i < evals1.size(); ++i)
    {
        phase1(i, i) = std::exp(-std::complex<double>(0, 1) * evals1(i));
    }

    Matrix Ustep1 = evecs1 * phase1 * evecs1.adjoint();

    Eigen::SelfAdjointEigenSolver<Matrix> eig2(K2);
    VectorReal evals2 = eig2.eigenvalues();
    Matrix evecs2 = eig2.eigenvectors();
    Matrix phase2 = Matrix::Zero(H.rows(), H.cols());
    for (int i = 0; i < evals2.size(); ++i)
    {
        phase2(i, i) = std::exp(-std::complex<double>(0, 1) * evals2(i));
    }

    Matrix Ustep2 = evecs2 * phase2 * evecs2.adjoint();


    U = Ustep1 * Ustep2 * U;
}

double spectralFunction(const VectorReal &eigVals, const Matrix &eigVecs, const Matrix &V, double omega, double sigma, double beta, double hbar)
{
    int D = eigVals.size();

    Matrix Vmn = eigVecs.adjoint() * V * eigVecs;

    double sum = 0.0;
    VectorReal boltzmannWeights(D);
    for(int n = 0; n < D; ++n)
    {
        boltzmannWeights(n) = std::exp(-eigVals(n) * beta);
    }

    boltzmannWeights /= boltzmannWeights.sum(); // Normalize the Boltzmann weights

    for (int n = 0; n < D; ++n)
    {
        for (int m = 0; m < D; ++m)
        {
            if (m == n) continue;
            double wmn = (eigVals(m) - eigVals(n))/hbar;

            if (std::abs(wmn) < TOL) continue;

            sum += std::norm(Vmn(m, n)) * std::exp(-(omega-wmn) * (omega-wmn)/(2.0 * sigma * sigma)) / (std::sqrt(2.0 * pi) * sigma) * (boltzmannWeights(n) + boltzmannWeights(m)) / 2.0;
        }
    }

    return sum;
}