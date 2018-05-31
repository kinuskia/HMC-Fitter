import numpy as np 
import matplotlib.pyplot as plt 

n, a0, a1, a2, a3, a4, a5, U, accept = np.loadtxt("data.txt", unpack = True)

keep = (U < 50)

parameters = [a0, a1, a2, a3, a4, a5]

plt.figure(1)
plt.hist(U[keep], 300)
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

counter_param = 0
counter_fig = 3
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















