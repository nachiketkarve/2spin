import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

S = 5
L = 1.0
kB = 1.0
Temp = 1.0
Ts = np.arange(0.01, 10.0, 0.01)

statistic = "expAvg"

TsPlot = np.array([])
Ws = np.array([])
vs = np.array([])


for T in Ts:
    FileName = "./../data/2spin_S" + str(S) + "_lambda" + "{:.6f}".format(L) + "_T" + "{:.6f}".format(T) + ".csv"
    try:
        df = pd.read_csv(FileName)
        DeltaE = df["deltaE"].values
        Einit = df["Ei"].values

        if statistic == "var":
            W = np.sum(DeltaE**2 * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp))) - (np.sum(DeltaE * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp))))**2
        elif statistic == "mean":
            W = np.sum(DeltaE * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp)))
        elif statistic == "expAvg":
            W = np.sum(np.exp(-DeltaE) * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp)))
        else:
            raise ValueError("Invalid statistic type. Choose from 'var', 'mean', or 'expAvg'.")

        #W = np.sum(DeltaE * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp)))

        #W = np.sum(np.exp(-DeltaE) * np.exp(-Einit / (kB * Temp))) / np.sum(np.exp(-Einit / (kB * Temp)))

        TsPlot = np.append(TsPlot, T)
        Ws = np.append(Ws, W)
        vs = np.append(vs, 2.0*L/(T*np.pi))

    except:
        continue

idx = np.argwhere(Ws > 0.0).flatten()
Ws = Ws[idx]
vs = vs[idx]
TsPlot = TsPlot[idx]
print(Ws)


mus = 1.0 / TsPlot
p = np.polyfit(np.log(mus[:]), np.log(Ws[:]/vs[:]**2), 1)

plt.figure()
#plt.scatter(mus, Ws,color='b')
plt.scatter(mus, Ws/vs**2,color='b')
if statistic == "var":
    plt.plot(mus, np.exp(p[1]) * mus**p[0], 'b--', label=r"$\frac{\text{Var}(\Delta K_0)}{\bar{v}^2} \sim \mu^{%.3f}$"%p[0])
elif statistic == "mean":
    plt.plot(mus, np.exp(p[1]) * mus**p[0], 'b--', label=r"$\frac{\langle \Delta K_0\rangle}{\bar{v}^2} \sim \mu^{%.3f}$"%p[0])
elif statistic == "expAvg":
    plt.plot(mus, np.exp(p[1]) * mus**p[0], 'b--', label=r"$\frac{\langle e^{-\Delta K_0}\rangle}{\bar{v}^2} \sim \mu^{%.3f}$"%p[0])

plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.show()