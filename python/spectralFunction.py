import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

S = 20
beta = 1e-2

file_name = "./../data/2spectralFunction_S" + "{:d}".format(S) + "_B" + "{:.6f}".format(beta) + ".csv"
df = pd.read_csv(file_name)

omega = df["frequency"].dropna().values
spectralFunctions = df["spectralFunction"].dropna().values

spectralFunctions = spectralFunctions / np.trapezoid(spectralFunctions, omega)

fig = plt.figure()
plt.plot(omega, spectralFunctions, '-o')
plt.xscale("log")
plt.yscale("log")
plt.ylim([1e-15, spectralFunctions.max() * 10])
plt.show()