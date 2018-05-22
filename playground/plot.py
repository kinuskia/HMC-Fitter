import numpy as np 
import matplotlib.pyplot as plt 

def exponential(x, N, tau):
	return N*np.exp(-x/tau)

q = np.loadtxt("sequence.txt", unpack = True)
g = np.loadtxt("gammas.txt", unpack = True)
t = np.loadtxt("times.txt", unpack = True)

t1 = np.arange(0, len(q),1)
t2 = np.arange(0, len(g),1)
t4 = np.arange(0, len(t),1)

plt.figure(1)
plt.plot(t1, q)
plt.savefig("sequence.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(2)
plt.plot(t2, g)
from scipy.optimize import curve_fit
popt, pcov = curve_fit(exponential, t2, g, p0 = [1, 10])
plt.plot(t2, exponential(t2, *popt))
#print("tau: ", popt[1], " +- ", np.sqrt(pcov[1][1]))
print("tau_int: ", -1./2+1./(1.-np.exp(-1/popt[1])))
plt.xlim(0, 20*popt[1])
plt.savefig("gammas.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(3)
plt.hist(q, 50)
plt.savefig("hist_sequence.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(4)
plt.plot(t4, t)
plt.xlim(0, 20*popt[1])
#print("W_theo = ", popt[1]/2*np.log(len(q)/popt[1]))
plt.savefig("times.pdf", format = "pdf", bbox_inches = "tight")
