import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#Import data
Data_File_Name = "data.dat"
colnames =  ["x", "R", "f'"] #name of columns
df = pd.read_csv("data.dat", sep="\t", header=None, names=colnames )
zeros = pd.read_csv("zeros.dat", sep="\t" )

#plotting  stuff
title = "Function Plot"
Plots_Files_Names = ["R_x1.pdf", "f_lambda.pdf"]
Plots_Titles = [r"Bessel Equation Solution with $\lambda=1$",
                r"Plot of $R(r=1; \lambda)=f(\lambda)$"]
Plots_Ylabel = [r"$R(x; \lambda=1)$", r"$f(\lambda)$"]
Plots_Xlabel = [r"$x$", r"\lambda"]
Plots_Xrange = [(0, 10), (0,15)]

print(zeros["a"])
for k in range(len(Plots_Files_Names)):
    plt.plot( df[colnames[0]], df[colnames[k+1]] )
    #Make x axis bold
    xlimit1, xlimit2= 0, 15 #x positions for the bold axis
    ylimit1, ylimit2 = 0, 0 #y positions for the bold axis
    #set the interval of the zeros
    zeros_midpoints = (zeros["a"]+ zeros["b"])/2
    n = len(zeros_midpoints)
    plt.scatter(zeros_midpoints, np.zeros(n), color="red", s=50,
            label= "Interval in which the zero lives", marker='o')
    plt.plot([xlimit1, xlimit2],
             [ylimit1, ylimit2], color="k", linewidth=1)
    plt.title(Plots_Titles[k])
    plt.xlabel( Plots_Xlabel[0] )
    plt.xlim(Plots_Xrange[k])
    plt.ylabel(Plots_Ylabel[k])
    plt.grid()
    plt.legend()
    plt.savefig(Plots_Files_Names[k])
    plt.close()
