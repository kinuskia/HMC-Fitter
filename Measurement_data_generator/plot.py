import numpy as np 
import matplotlib.pyplot as plt 

x, y, dy = np.loadtxt("measurements.txt", unpack = True)

def fit (x, a, b, c, d, e):
	return a + b * np.exp(-c*x)*np.sin(d*x+e) 
from scipy.optimize import curve_fit
popt, pcov = curve_fit(fit, x, y, p0= [2.3, 1.2, 0.3, 2.0, 0.8], sigma=dy, absolute_sigma = True)
print(popt[0], " + - ", np.sqrt(pcov[0][0]))
print(popt[1], " + - ", np.sqrt(pcov[1][1]))
print(popt[2], " + - ", np.sqrt(pcov[2][2]))
print(popt[3], " + - ", np.sqrt(pcov[3][3]))
print(popt[4], " + - ", np.sqrt(pcov[4][4]))



plt.errorbar(x, y, yerr = dy, linestyle = "none", capsize = 2)
values = np.linspace(0, 10, 300)
plt.plot(values, fit(values, *popt))
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("measurements.pdf", format = "pdf", bbox_inches = "tight")

chi2 = np.sum((y-fit(x, *popt))**2/dy/dy)
d_of_freedom = len(x) - len(popt)
print(chi2/d_of_freedom)

