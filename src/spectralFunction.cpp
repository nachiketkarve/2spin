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
    double beta = 1.0;

    double sigma = 1.0;

    double freqMin = 0.001;
    double freqMax = 10.0;
    int freqNum = 1000;

    if (argc != 4)
    {
        throw std::invalid_argument("usage: spectralFunction <freqMin> <freqMax> <freqNum>");
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
    beta = params["beta"];
    sigma = params["sigma"];

    freqMin = std::atof(argv[1]);
    freqMax = std::atof(argv[2]);
    freqNum = std::atoi(argv[3]);

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
    if (freqMin <= 0.0 || freqMax <= 0.0 || freqNum <= 0 || freqMin >= freqMax)
    {
        throw std::invalid_argument("frequency parameters not valid. Input [freqMin], [freqMax], [freqNum] with 0 < freqMin < freqMax and freqNum > 0.");
    }
    if (beta < 0.0)
    {
        throw std::invalid_argument("beta must be non-negative");
    }

    Jx = J[0];
    Jy = J[1];
    Jz = J[2];
    Ax = A[0];
    Ay = A[1];
    Az = A[2];

    std::string fileName = saveFolder + "2spectralFunction_S" + std::to_string(S) + "_B" + std::to_string(beta) + ".csv";

    int dim = (2 * S + 1) * (2 * S + 1);
    Matrix H(dim, dim);
    Matrix S1zS2z(dim, dim);
    constructHamiltonian(H, S, S, Jx, Jy, Jz, Ax, Ay, Az, hbar);
    constructS1zS2z(S1zS2z, S, S, hbar);

    Eigen::SelfAdjointEigenSolver<Matrix> eigSolver(H);
    VectorReal eigVals = eigSolver.eigenvalues();
    Matrix eigVecs = eigSolver.eigenvectors();

    VectorReal freqs(freqNum);
    VectorReal exponents = VectorReal::LinSpaced(freqNum, std::log(freqMin), std::log(freqMax));
    VectorReal spectralFunctions(freqNum);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < freqNum; ++i)
    {
        freqs(i) = std::exp(exponents(i));
        spectralFunctions(i) = spectralFunction(eigVals, eigVecs, S1zS2z, freqs(i), sigma, beta, hbar);
    }

    std::ofstream file;
    file.open(fileName);
    file << std::setprecision(15);
    file << "frequency" << "," << "spectralFunction";
    for (int i = 0; i < freqNum; ++i)
    {
        file << "\n" << freqs(i) << "," << spectralFunctions(i);
    }

    file.close();

}