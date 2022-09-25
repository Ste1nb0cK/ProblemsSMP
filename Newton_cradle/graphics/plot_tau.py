import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data1 = pd.read_table("data_1.txt", header = None)
data2 = pd.read_table("data_2.txt", header = None)
data3 = pd.read_table("data_3.txt", header = None)
data4 = pd.read_table("data_4.txt", header = None)
data5 = pd.read_table("data_5.txt", header = None)
data6 = pd.read_table("data_6.txt", header = None)
data7 = pd.read_table("data_7.txt", header = None)

"""plt.plot(data1[0], data1[1], label = "$K = 0.1e10$")
plt.plot(data2[0], data2[1], label = "$K = 0.2e10$")
plt.plot(data3[0], data3[1], label = "$K = 0.5e10$")
plt.plot(data4[0], data4[1], label = "$K = 1e10$")
plt.plot(data5[0], data5[1], label = "$K = 2e10$")
plt.plot(data6[0], data6[1], label = "$K = 5e10$")
plt.plot(data7[0], data7[1], label = "$K = 10e10$")
plt.title("$\\tau$ vs $t$")
plt.ylabel("$\\tau$ (N m)")
plt.xlabel("t (s)")
plt.xlim(0.173, 0.178)
plt.legend()
plt.grid()
plt.savefig("tau_vs_t.pdf")"""

t0 = 0.17457
beta = -0.403
alpha = 0.400
data6[0]=5e10**(-beta)*(data6[0]-t0)
data6[1]=5e10**(-alpha)*data6[1]

plt.plot(0.1e10**(-beta)*(data1[0]-t0), 0.1e10**(-alpha)*data1[1], label = "$K = 0.1e10$")
plt.plot(0.2e10**(-beta)*(data2[0]-t0), 0.2e10**(-alpha)*data2[1], label = "$K = 0.2e10$")
plt.plot(0.5e10**(-beta)*(data3[0]-t0), 0.5e10**(-alpha)*data3[1], label = "$K = 0.5e10$")
plt.plot(1e10**(-beta)*(data4[0]-t0), 1e10**(-alpha)*data4[1], label = "$K = 1e10$")
plt.plot(2e10**(-beta)*(data5[0]-t0), 2e10**(-alpha)*data5[1], label = "$K = 2e10$")
plt.plot(data6[0], data6[1], label="$K = 1e10$")
plt.plot(10e10**(-beta)*(data7[0]-t0),10e10**(-alpha)*data7[1], label = "$K = 10e10$")
plt.title("$K^{-\\alpha}\\tau$ vs $K^{-\\beta}(t-t_{0})$")
plt.ylabel("$K^{-\\alpha}\\tau$")
plt.xlabel("$K^{-\\beta}(t-t_{0})$")
#plt.xlim(0.173, 0.178)
plt.xlim(0, 20)
plt.legend()
plt.grid()
plt.savefig("colapso.pdf")
