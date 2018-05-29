import numpy as np 
import matplotlib.pyplot as plt 

x, y, dy = np.loadtxt("correlator_Ps-Ps", unpack = True)


plt.errorbar(x, y, yerr = dy, linestyle = "none", capsize = 2)

plt.savefig("measurements.pdf", format = "pdf", bbox_inches = "tight")


