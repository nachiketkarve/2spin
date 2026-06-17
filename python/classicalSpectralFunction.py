import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

lmd = 0.0
T = 1000.0
beta = 0.01


file_name = "./../data/classicalSpectralFunction_lambda" + "{:.6f}".format(lmd) + "_T" + "{:.6f}".format(T) + "_B" + "{:.6f}".format(beta) + ".csv"
df = pd.read_csv(file_name)

frequency = df["frequency"].dropna().values
SpectralFunction = df["SpectralFunction"].dropna().values

SpectralFunction = SpectralFunction / np.trapezoid(SpectralFunction, frequency)

fitIdx = np.arange(1, len(frequency)*0.001, 1, dtype=int)
p = np.polyfit(np.log(frequency[fitIdx]), np.log(SpectralFunction[fitIdx]), 1)

fig = plt.figure()
plt.plot(frequency[1:], SpectralFunction[1:], '-o')
plt.plot(frequency[fitIdx], np.exp(p[1]) * frequency[fitIdx] ** p[0], '--', label=r"$S(\omega)\sim \omega^{%.3f}$"%p[0])
plt.xlabel("Frequency")
plt.ylabel("Spectral Function")
plt.yscale("log")
plt.xscale("log")
plt.legend()

plt.show()