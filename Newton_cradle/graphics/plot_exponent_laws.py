import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
plt.figure()
K = 1e10*np.array([0.1, 0.2, 0.5, 1, 2, 5, 10])
tau_max = 1e7*np.array([2.85, 3.75, 5.41, 7.11, 9.445, 13.6, 18.02])
K = np.log10(K)
tau_max = np.log10(tau_max)
DT = 1e-3*np.array([1.55, 1.15, 0.8, 0.6, 0.45, 0.3, 0.25])
DT = np.log10(DT)
plt.plot(K, tau_max)
plt.scatter(K, tau_max)
plt.title("$\\tau_{max}$ vs K")
plt.ylabel("$log_{10}(\\tau_{max})$")
plt.xlabel("$log_{10}(K)$")
plt.grid()
plt.savefig("tau_max.pdf")

plt.figure()

plt.plot(K, DT)
plt.scatter(K, DT)
plt.title("$t_{max}$ vs K")
plt.ylabel("$log_{10}(t_{max})$")
plt.xlabel("$log_{10}(K)$")
plt.grid()
plt.savefig("t_max.pdf")
