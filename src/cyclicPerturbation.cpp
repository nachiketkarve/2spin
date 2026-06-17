#include "2spin.hpp"
#include "nlohmann/json.hpp"
#include <cstdlib>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int S = 4;
    double hbar = 1.0;
    double Jx = 1.0, Jy = 1.0, Jz = 0.5;
    double Ax = 0.0, Ay = 0.0, Az = 0.0;
    double dt = 0.1;
    double lambda = 0.1;
    double T = 1.0;

    if (argc != 3)
    {
        throw std::invalid_argument("usage: 2spin <maxTime> <lambda>");
    }

    std::ifstream dataFile("_params.json");
    if (!dataFile)
    {
        throw std::runtime_error("could not open _params.json");
    }

    json params = json::parse(dataFile);
    S = params["S"];
    hbar = params["hbar"];
    std::vector<double> J = params["J"];
    std::vector<double> A = params["A"];
    std::string saveFolder = params["saveFolder"];
    int integratorOrder = params["integratorOrder"];
    dt = params["dt"];

    T = std::atof(argv[1]);
    lambda = std::atof(argv[2]);
    if (S < 0)
    {
        throw std::invalid_argument("S must be non-negative");
    }
    if (J.size() != 3 || A.size() != 3)
    {
        throw std::invalid_argument("J and A must each contain exactly 3 values");
    }
    if (hbar == 0.0)
    {
        throw std::invalid_argument("hbar must be non-zero");
    }
    if (T <= 0.0 || dt <= 0.0)
    {
        throw std::invalid_argument("T and dt must be positive");
    }

    Jx = J[0];
    Jy = J[1];
    Jz = J[2];
    Ax = A[0];
    Ay = A[1];
    Az = A[2];

    int Nsteps = static_cast<int>(T / dt);
    if (Nsteps <= 0)
    {
        throw std::invalid_argument("T must be at least one dt");
    }
    dt = T / Nsteps; // Adjust dt to fit exactly into T
    double mu = 2.0 * pi / T;

    std::string fileName = saveFolder + "2spin_S" + std::to_string(S) + "_lambda" + std::to_string(lambda) + "_T" + std::to_string(T) + ".csv";

    int dim = (2 * S + 1) * (2 * S + 1);
    Matrix H(dim, dim);
    Matrix S1zS2z(dim, dim);
    constructHamiltonian(H, S, S, Jx, Jy, Jz, Ax, Ay, Az, hbar);
    constructS1zS2z(S1zS2z, S, S, hbar);

    Eigen::SelfAdjointEigenSolver<Matrix> eigSolver(H);
    VectorReal eigVals = eigSolver.eigenvalues();
    Matrix eigVecs = eigSolver.eigenvectors();

    Matrix U = Matrix::Identity(dim, dim);
    for (int step = 0; step < Nsteps; step++)
    {
        if (integratorOrder == 4)
        {
            UpdateUnitary4thOrder(H, S1zS2z, U, hbar, lambda, mu, step * dt, dt);
        }
        else
        {
            UpdateUnitary(H, S1zS2z, U, hbar, lambda, mu, step * dt, dt);
        }
    }

    Matrix I = Matrix::Identity(dim, dim);
    std::cout << "Unitarity error: " << (U.adjoint() * U - I).norm() << std::endl;

    Matrix H_eig = eigVecs.adjoint() * H * eigVecs;
    
    Matrix D = Matrix::Zero(dim, dim);
    for (int i = 0; i < dim; i++)
    {
        D(i, i) = std::complex<double>(eigVals(i), 0.0);
    }

    double diagonalization_error = (H_eig - D).norm();
    std::cout << "diagonalization_error = " << diagonalization_error << std::endl;

    VectorReal Es(dim);
    for (int i = 0; i < dim; i++)
    {
        
        Vector psi_i = eigVecs.col(i);
        Vector psi_f = U * psi_i;

        double norm = std::real(psi_f.dot(psi_f));
        double Ef = std::real(psi_f.dot(H * psi_f)) / norm;

        Es(i) = Ef;
    }

    Vector psi0 = eigVecs.col(0);
    Vector psi0f = U * psi0;

    double norm0 = std::real(psi0f.dot(psi0f));
    double Ef0_direct = std::real(psi0f.dot(H * psi0f)) / norm0;
    double ground_delta_direct = Ef0_direct - eigVals(0);

    std::cout << std::setprecision(17);

    std::cout << "E0 = " << eigVals(0) << std::endl;
    std::cout << "Ef0 direct = " << Ef0_direct << std::endl;
    std::cout << "Ef0 - E0 direct = " << ground_delta_direct << std::endl;
    std::cout << "norm0 direct = " << norm0 << std::endl;


    std::ofstream file;
    file.open(fileName);
    if (!file)
    {
        throw std::runtime_error("could not open output file: " + fileName);
    }
    file << std::setprecision(15);

    file << "Ei,Ef,deltaE";
    for (int i = 0; i < dim; i++)
    {
        file << "\n" << eigVals(i) << "," << Es(i) << "," << Es(i) - eigVals(i);
    }

    file.close();





}
