import numpy as np 
import matplotlib.pyplot as plt 

n, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, U, accept = np.loadtxt("data.txt", unpack = True)
#a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, U, accept = np.loadtxt("LPT-Cluster/data1.txt", unpack = True)
#n = np.arange(0, len(a0))
keep = (U <100)

parameters = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25]

T = 1e-1
p = 26

plt.figure(1)
plt.hist(U[keep], 300)
from scipy.optimize import curve_fit
plt.xlabel("$\\chi^2_\\mathrm{red}$")
plt.ylabel("#")
plt.savefig("plots/hist_U.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)
plt.figure(2)
plt.plot(n[keep], U[keep])
plt.xlabel("$n$")
plt.ylabel("$\\chi^2_\\mathrm{red}$")
plt.savefig("plots/U.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(3)
plt.plot(n[keep], accept[keep])
plt.xlabel("$n$")
plt.ylabel("acceptance rate")
plt.savefig("plots/acceptance.pdf", format = "pdf", bbox_inches = "tight")
plt.close(3)

counter_param = 0
counter_fig = 4
n = n[keep]

for item in parameters:
	item = item[keep]
	plt.figure(counter_fig)
	plt.plot(n, item)
	plt.xlabel("n")
	lab1 = "popt[" + str(counter_param) + "]"
	plt.ylabel(lab1)
	filename1 = "plots/popt" + str(counter_param) + ".pdf"
	plt.savefig(filename1, format = "pdf", bbox_inches = "tight")
	plt.close(counter_fig)
	counter_fig = counter_fig + 1
	
	plt.figure(counter_fig)
	plt.hist(item, 300)
	plt.xlabel(lab1)
	plt.ylabel("#")
	filename2 = "plots/hist_popt" + str(counter_param) + ".pdf"
	plt.savefig(filename2, format = "pdf", bbox_inches = "tight")
	plt.close(counter_fig)
	counter_fig = counter_fig + 1
	counter_param = counter_param + 1

mean = np.mean(U[keep])
chi2redmin = np.min(U[keep])
diff_theo = p*T/2
print((mean-chi2redmin)/diff_theo )















