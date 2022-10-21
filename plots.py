import numpy as np
import matplotlib.pyplot as plt
import functions
from functions import *

# -------------data--------------------

m = np.array([[5,1], [3,2]])

tau_avg = np.genfromtxt("Tau_avg_5_2_3_1")
tau_avg_A = np.genfromtxt("Tau_avg_A_5_2_3_1")
theory_tau = np.genfromtxt("Theory_tau_5_1_3_2.txt")

beta = np.linspace(0, 0.01 ,11)
beta2 = np.linspace(0, 0.012, 13)

N=100

# -------------------------------------


fig, axs = plt.subplots(2, sharex=True)
fig.suptitle('Fixation Times in Evolutionary Games under Weak Selection')

axs[0].plot(beta, tau_avg, linestyle='', marker='D', markerfacecolor="None",  markeredgecolor="Black")
axs[0].plot(beta2, get_theory_tau(N, beta2, m), color="Black")
axs[0].plot(beta2, weak_fixation_time(N, m, beta2), linestyle="dashed", color='cyan')

axs[0].set_ylim(0.895, 1.005)  
axs[0].set_xlim(-0.0005, 0.0105)
axs[0].set_yticks(np.arange(0.9, 1, 0.025), minor=True)
axs[0].tick_params(right=True, top=True)
axs[0].set_ylabel(r"$\tau_{1}$")

axs[1].plot(beta, tau_avg_A, linestyle='', marker='o', markerfacecolor="None",  markeredgecolor="Black")
axs[1].plot(beta2, theory_tau, color="Black")
axs[1].plot(beta2, weak_fixation_time_A(N, beta2, m), linestyle="dashed", color='cyan')

axs[1].set_ylim(0.895, 1.005)
axs[1].set_xlim(-0.0005, 0.0105)
axs[1].set_yticks(np.arange(0.9, 1, 0.025), minor=True)
axs[1].set_xticks(np.arange(0, 0.01, 0.001), minor=True)
axs[1].tick_params(right=True, top=True)

axs[1].set_xlabel(r"$\beta$")    
axs[1].set_ylabel(r"$\tau^{A}_{1}$")

plt.show()

