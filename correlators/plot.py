import numpy as np 
import matplotlib.pyplot as plt 

x, y, dy = np.loadtxt("correlator_Pl-Pl", unpack = True)

T = 64.
def fit(t, m1, m2, Z1, Z2, Ar, mr):
	result = Z1*Z1/2./m1*(np.exp(-m1*t)+np.exp(-m1*(T-t)))
	result += Z2*Z2/2./m2*(np.exp(-m2*t)+np.exp(-m2*(T-t)))
	result += Ar*(np.exp(-mr*t)+np.exp(-mr*(T-t)))
	return result

from scipy.optimize import curve_fit
for i in range(0, len(x)):
	print(x[i], y[i], dy[i])

popt, pcov = curve_fit(fit, x, y, p0 = [2.9, 3.6, 8, 1.4, 0.05, 1.05], sigma = dy, absolute_sigma=True) 
#p0 = [3.59, 1.638, 2.4358, 1.028, 0.04083, 0.9912]
plt.errorbar(x, y, yerr = dy, linestyle = "none", capsize = 2)
values = np.linspace(x[0], x[-1],100)
plt.yscale("log")
plt.plot(fit(values, *popt))

plt.savefig("measurements.pdf", format = "pdf", bbox_inches = "tight")
chi2 = np.sum((y-fit(x,*popt))**2/dy**2)
dof = len(x) - len(popt)
for i in range(0, len(popt)):
	print(i, popt[i], "+-", np.sqrt(pcov[i][i]))

print("chi2/dof: ", chi2/dof)