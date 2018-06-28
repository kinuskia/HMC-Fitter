import numpy as np 
import matplotlib.pyplot as plt 
import os

# Data is collected in these arrays
a0 = np.empty((0))
a1 = np.empty((0))
a2 = np.empty((0))
a3 = np.empty((0))
a4 = np.empty((0))
a5 = np.empty((0))
a6 = np.empty((0))
a7 = np.empty((0))
a8 = np.empty((0))
a9 = np.empty((0))
a10 = np.empty((0))
a11 = np.empty((0))
a12 = np.empty((0))
a13 = np.empty((0))
a14 = np.empty((0))
a15 = np.empty((0))
a16 = np.empty((0))
a17 = np.empty((0))
a18 = np.empty((0))
a19 = np.empty((0))
a20 = np.empty((0))
a21 = np.empty((0))
a22 = np.empty((0))
a23 = np.empty((0)) 
a24 = np.empty((0)) 
a25 = np.empty((0)) 
a26 = np.empty((0)) 
a27 = np.empty((0)) 
a28 = np.empty((0)) 
a29 = np.empty((0)) 
U = np.empty((0))
accept = np.empty((0))


# Fill arrays with data
counter_files = 0
for file in os.listdir("/Users/Kianusch/Documents/Studium/Studiensemester/SoSe18/stage/HMC-Fitter/LPT-Cluster"):
	filename = os.fsdecode(file)
	if filename.endswith(".txt"):
		counter_files = counter_files + 1
		b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, bU, baccept = np.loadtxt("LPT-Cluster/" + filename, unpack = True)
		a0 = np.append(a0, b0)
		a1 = np.append(a1, b1)
		a2 = np.append(a2, b2)
		a3 = np.append(a3, b3)
		a4 = np.append(a4, b4)
		a5 = np.append(a5, b5)
		a6 = np.append(a6, b6)
		a7 = np.append(a7, b7)
		a8 = np.append(a8, b8)
		a9 = np.append(a9, b9)
		a10 = np.append(a10, b10)
		a11 = np.append(a11, b11)
		a12 = np.append(a12, b12)
		a13 = np.append(a13, b13)
		a14 = np.append(a14, b14)
		a15 = np.append(a15, b15)
		a16 = np.append(a16, b16)
		a17 = np.append(a17, b17)
		a18 = np.append(a18, b18)
		a19 = np.append(a19, b19)
		a20 = np.append(a20, b20)
		a21 = np.append(a21, b21)
		a22 = np.append(a22, b22)
		a23 = np.append(a23, b23)
		a24 = np.append(a24, b24)
		a25 = np.append(a25, b25)
		a26 = np.append(a26, b26)
		a27 = np.append(a27, b27)
		a28 = np.append(a28, b28)
		a29 = np.append(a29, b29)
		U = np.append(U, bU)
		accept = np.append(accept, baccept)
n = np.arange(0, len(a0), 1)


keep = (U <42)

parameters = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29]

T = 5e-3
p = 30.

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

mean = np.mean(U[keep])
chi2redmin = np.min(U[keep])
diff_theo = p*T/2
print((mean-chi2redmin)/diff_theo )















