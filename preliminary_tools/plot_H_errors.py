import numpy as np 
import matplotlib.pyplot as plt 

H, error = np.loadtxt("H_errors.txt", unpack = True)


plt.figure(1)
plt.scatter(H, abs(error))
plt.xlabel("H")
plt.ylabel("relative discretization_error")
plt.xscale("log")
plt.yscale("log")
plt.savefig("discretization_error.pdf", format = "pdf", bbox_inches = "tight")





