import numpy as np 
import matplotlib.pyplot as plt 

x, y, dy = np.loadtxt("correlator_Pl-Pl", unpack = True, skiprows = 1)
x2, y2, dy2 = np.loadtxt("correlator_Ps-A0l", unpack = True, skiprows = 1)
x3, y3, dy3 = np.loadtxt("correlator_A0l-Ps", unpack = True, skiprows = 1)

# T = 64.
# def fit(t, m1, m2, Z1, Z2, Ar, mr):
# 	result = Z1*Z1/2./m1*(np.exp(-m1*t)+np.exp(-m1*(T-t)))
# 	result += Z2*Z2/2./m2*(np.exp(-m2*t)+np.exp(-m2*(T-t)))
# 	result += Ar*(np.exp(-mr*t)+np.exp(-mr*(T-t)))
# 	return result

from scipy.optimize import curve_fit

#popt, pcov = curve_fit(fit, x, y, p0 = [0.99, 1.64, 0.28, 1.03, 0.83, 3.58], sigma = dy, absolute_sigma=True) 
#p0 = [0.99, 1.64, 0.28, 1.03, 0.83, 3.58]
#plt.errorbar(x, y, yerr = dy, linestyle = "none", capsize = 2, label="Pl-Pl")
#plt.errorbar(x2, y2, yerr = dy2, linestyle = "none", capsize = 2, label="A0l-A0s")
#plt.errorbar(x3, y3, yerr = dy3, linestyle = "none", capsize = 2, label="A0s-A0l")
plt.errorbar(x3, (y3-y2)/y2, yerr = np.sqrt(dy3*dy3+dy2*dy2)/y2, linestyle = "none", capsize = 2, label="difference")
values = np.linspace(0., 29., num = 100)
#plt.yscale("log")
plt.legend(loc = "best")
#plt.plot(values, fit(values, *popt))

plt.savefig("measurements.pdf", format = "pdf", bbox_inches = "tight")
#chi2 = np.sum((y-fit(x,*popt))**2/dy**2)
#dof = len(x) - len(popt)
#for i in range(0, len(popt)):
#	print(i, popt[i], "+-", np.sqrt(pcov[i][i]))

#print("chi2/dof: ", chi2/dof)