#include "classical2spin.hpp"
#include "nlohmann/json.hpp"
#include <random>
#include <unsupported/Eigen/FFT>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    double Jx = 1.0, Jy = 1.0, Jz = 0.5;
    double Ax = 0.0, Ay = 0.0, Az = 0.0;
    double dt = 0.1;
    double lambda = 0.1;
    double T = 1.0;
    double beta = 1.0;
    int realizations = 10;

    if (argc != 4)
    {
        throw std::invalid_argument("usage: 2spin <maxTime> <lambda> <seed>");
    }

    std::ifstream dataFile("_params.json");
    json params = json::parse(dataFile);
    std::vector<double> J = params["J"];
    std::vector<double> A = params["A"];
    realizations = params["realizations"];
    std::string saveFolder = params["saveFolder"];
    dt = params["dt"];
    beta = params["beta"];

    T = std::atof(argv[1]);
    lambda = std::atof(argv[2]);
    int seed = std::atoi(argv[3]);

    Jx = J[0];
    Jy = J[1];
    Jz = J[2];
    Ax = A[0];
    Ay = A[1];
    Az = A[2];

    int Nsteps = static_cast<int>(T / dt);
    dt = T / Nsteps; // Adjust dt to fit exactly into T
    double mu = 2.0 * pi / T;

    std::string fileName = saveFolder + "classicalSpectralFunction_lambda" + std::to_string(lambda) + "_T" + std::to_string(T) + "_B" + std::to_string(beta) + ".csv";


    std::mt19937                        generator(seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::uniform_real_distribution<double> uniformPhi(0.0, 2.0 * pi);

    Eigen::FFT<double> fft;

    VectorReal frequencies = VectorReal::LinSpaced(Nsteps, 0.0, Nsteps - 1.0) / (Nsteps * dt);
    VectorReal SpectralFnAvg(Nsteps);
    SpectralFnAvg.setZero();
    double Z = 0.0;

    for (int realization = 0; realization < realizations; realization++)
    {

        double z1 = 2.0 * uniform01(generator) - 1.0;
        double phi1 = uniformPhi(generator);
        double r1 = std::sqrt(1.0 - z1 * z1);

        double z2 = 2.0 * uniform01(generator) - 1.0;
        double phi2 = uniformPhi(generator);
        double r2 = std::sqrt(1.0 - z2 * z2);


        VectorReal S1(3), S2(3);
        S1(0) = r1 * std::cos(phi1);
        S1(1) = r1 * std::sin(phi1);
        S1(2) = z1;
        S2(0) = r2 * std::cos(phi2);
        S2(1) = r2 * std::sin(phi2);
        S2(2) = z2;

        VectorReal Observable(Nsteps);
        Vector currSpectrum(Nsteps);

        for (int i = 0; i < Nsteps; i++)
        {
            double t = i * dt;
            Observable(i) = S1(2) * S2(2);
            Update(S1, S2, Jx, Jy, Jz, Ax, Ay, Az, lambda, mu, t, dt);
        }

        Observable.array() -= Observable.mean();

        fft.fwd(currSpectrum, Observable);
        SpectralFnAvg += currSpectrum.cwiseAbs2()/realizations * std::exp(-beta * Hamiltonian(S1, S2, Jx, Jy, Jz, Ax, Ay, Az, lambda, mu, T));
        Z += std::exp(-beta * Hamiltonian(S1, S2, Jx, Jy, Jz, Ax, Ay, Az, lambda, mu, T));
    }

    SpectralFnAvg /= Z;

    std::ofstream file;
    file.open(fileName);
    file << std::setprecision(15);

    file << "frequency,SpectralFunction";
    for (int i = 0; i < int(Nsteps/2); i++)
    {
        file << "\n" << frequencies(i) << "," << SpectralFnAvg(i);
    }

    file.close();





}