import numpy as np 
import matplotlib.pyplot as plt 

n_steps = np.loadtxt("correlation_times.txt", unpack = True)



plt.figure(1)
plt.hist(n_steps, 50)
plt.xlabel("correlation time")
plt.ylabel("#")
plt.savefig("corr_lengths.pdf", format = "pdf", bbox_inches = "tight")





