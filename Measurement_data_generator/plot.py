import numpy as np 
import matplotlib.pyplot as plt 

x, y, dy = np.loadtxt("measurements.txt", unpack = True)

def fit (x, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9):
	return a0 * np.arctan(a1*x+a2) + a3 + a4*np.exp(-(x-a5)*(x-a5)/2./a6/a6) +a7*np.exp(-(x-a8)*(x-a8)/2./a9/a9) 
from scipy.optimize import curve_fit
popt, pcov = curve_fit(fit, x, y, p0= [-1.0, 1.2, -3.0, 3.0, 4., 4., 1.2, 5., 7., 0.7], sigma=dy, absolute_sigma = True)
print(popt[0], " + - ", np.sqrt(pcov[0][0]))
print(popt[1], " + - ", np.sqrt(pcov[1][1]))
print(popt[2], " + - ", np.sqrt(pcov[2][2]))
print(popt[3], " + - ", np.sqrt(pcov[3][3]))
print(popt[4], " + - ", np.sqrt(pcov[4][4]))
print(popt[5], " + - ", np.sqrt(pcov[5][5]))
print(popt[6], " + - ", np.sqrt(pcov[6][6]))
print(popt[7], " + - ", np.sqrt(pcov[7][7]))
print(popt[8], " + - ", np.sqrt(pcov[8][8]))
print(popt[9], " + - ", np.sqrt(pcov[9][9]))



plt.errorbar(x, y, yerr = dy, linestyle = "none", capsize = 2)
values = np.linspace(0, 10, 300)
plt.plot(values, fit(values, *popt))
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("measurements.pdf", format = "pdf", bbox_inches = "tight")

chi2 = np.sum((y-fit(x, *popt))**2/dy/dy)
d_of_freedom = len(x) - len(popt)
print(chi2/d_of_freedom)

