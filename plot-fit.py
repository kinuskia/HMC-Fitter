import numpy as np 
import matplotlib.pyplot as plt 

n, a, b, c, d, U, accept = np.loadtxt("data.txt", unpack = True)

plt.figure(1)
plt.plot(n, a)
plt.xlabel("n")
plt.ylabel("popt[0]")
plt.savefig("plots/popt0.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(2)
plt.plot(n, b)
plt.xlabel("n")
plt.ylabel("popt[1]")
plt.savefig("plots/popt1.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(3)
plt.plot(n, c)
plt.xlabel("n")
plt.ylabel("popt[2]")
plt.savefig("plots/popt2.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(4)
plt.plot(n, d)
plt.xlabel("n")
plt.ylabel("popt[3]")
plt.savefig("plots/popt3.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(5)
plt.hist(a, 1000)
plt.xlabel("popt[0]")
plt.ylabel("#")
plt.savefig("plots/hist_popt0.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(6)
plt.hist(b, 1000)
plt.xlabel("popt[1]")
plt.ylabel("#")
plt.savefig("plots/hist_popt1.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(7)
plt.hist(c, 1000)
plt.xlabel("popt[2]")
plt.ylabel("#")
plt.savefig("plots/hist_popt2.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(8)
plt.hist(d, 1000)
plt.xlabel("popt[3]")
plt.ylabel("#")
plt.savefig("plots/hist_popt3.pdf", format = "pdf", bbox_inches = "tight")

plt.figure(9)
plt.hist(U, 1000)
plt.xlabel("$\\chi^2_\\mathrm{red}$")
plt.ylabel("#")
plt.savefig("plots/hist_U.pdf", format = "pdf", bbox_inches = "tight")



