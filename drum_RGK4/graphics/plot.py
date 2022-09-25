import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data1 = pd.read_table("data_1.txt", header = None)
data2 = pd.read_table("data_2.txt", header = None)
data3 = pd.read_table("data_3.txt", header = None)
data4 = pd.read_table("data_4.txt", header = None)
data5 = pd.read_table("data_5.txt", header = None)

plt.plot(data1[0], data1[1], label = "$\lambda = 2.43797$")
plt.plot(data2[0], data2[1], label = "$\lambda = 5.4858$")
plt.plot(data3[0], data3[1], label = "$\lambda = 8.66974$")
plt.plot(data4[0], data4[1], label = "$\lambda = 11.7982$")
plt.plot(data5[0], data5[1], label = "$\lambda = 14.9112$")
plt.title("Theoritical Normal modes of vibration")
plt.xlabel("r")
plt.ylabel("R(r)")
plt.xlim(0, 1)
plt.legend()
plt.grid()
plt.savefig("Normal_modes_theo.pdf")
