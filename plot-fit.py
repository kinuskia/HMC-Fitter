import numpy as np 
import matplotlib.pyplot as plt 

n, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, U, accept = np.loadtxt("data.txt", unpack = True)
#a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, U, accept = np.loadtxt("LPT-Cluster/data1.txt", unpack = True)
#n = np.arange(0, len(a0))
keep = (U <40.85619)

parameters = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23]

T = 1e-6
p = 24

from scipy.special import factorial
def chi2_dist(z, chi2min, f_2, beta, N):
	if z > chi2min:
		return N*beta**(f_2)/factorial(f_2-1)*np.exp(-beta*(z-chi2min))*(z-chi2min)**(f_2-1)
	else:
		0
def gauss(z, mu, sigma, f, beta, N):
	sigma = sigma*np.sqrt(f/2./beta)
	return N/np.sqrt(2*np.pi*sigma*sigma)*np.exp(-(z-mu)**2/2/sigma/sigma)
n_bins = 128
f = 440
plt.figure(1)
n_chi2, bins_chi2, patches_chi2 = plt.hist(U[keep], n_bins, normed = True) 
xvalues_chi2 = np.linspace(40.856155, 40.856182, 200)
yvalues_chi2 = np.zeros(len(xvalues_chi2))
for i in range (0, len(yvalues_chi2)):
	yvalues_chi2[i] = chi2_dist(xvalues_chi2[i], 40.8561529, 12, 1e6, 1)
plt.plot(xvalues_chi2, yvalues_chi2)
plt.xlabel("$\\chi^2_\\mathrm{red}$")
plt.ylabel("density")
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
	n_popt, bins_popt, patches_popt = plt.hist(item, n_bins, normed = True)
	if counter_param == 0:
		xvalues_popt = np.linspace(0.985310, 0.985320, 200)
		yvalues_popt = np.zeros(len(xvalues_popt))
		yvalues_popt = gauss(xvalues_popt, 9.85315322e-1, 8.201e-5, f, 1e6, 1)
		plt.plot(xvalues_popt, yvalues_popt)
	if counter_param == 1:
		xvalues_popt = np.linspace(1.35160, 1.35180, 200)
		yvalues_popt = np.zeros(len(xvalues_popt))
		yvalues_popt = gauss(xvalues_popt, 1.35170252, 1.679e-3, f, 1e6, 1)
		plt.plot(xvalues_popt, yvalues_popt)
	plt.xlabel(lab1)
	plt.ylabel("density")
	filename2 = "plots/hist_popt" + str(counter_param) + ".pdf"
	plt.savefig(filename2, format = "pdf", bbox_inches = "tight")
	plt.close(counter_fig)
	counter_fig = counter_fig + 1
	counter_param = counter_param + 1

mean = np.mean(U[keep])
chi2redmin = np.min(U[keep])
diff_theo = p*T/2
print((mean-chi2redmin)/diff_theo )














