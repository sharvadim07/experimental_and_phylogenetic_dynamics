import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import *
from matplotlib.backends.backend_pdf import PdfPages

'''
	1st argument -- input filename
	2nd argument -- output filename
        3rd argument -- perform log-scaling (YES/NO)
	4th argument -- file with experimental points
'''

if sys.argv[3] == "YES":
	log_scale = 1
else:
	log_scale=0	

data = np.loadtxt(sys.argv[1])
if log_scale == 1:
	data [: , :] = data [: , :] + 1
	data [:, 0] = data [:, 0] - 1

try:
	exp_data_ = np.loadtxt(sys.argv[4], dtype=object)
	exp_times = list(map(float, exp_data_[1:, 0]))
	exp_data = {}
	for i, func in enumerate(exp_data_[0][1:]):
		col = exp_data_[1:, i+1]
		exp_data[func] = list(map(float, col))
		if log_scale == 1:
			exp_data[func] = [x+1 for x in exp_data[func]]

except:
	exp_data = {}


colors = [(0, 160, 0), (0, 0, 200), (200, 0, 0), (200, 200, 0), (255, 128, 0),
          (200, 0, 0), (200, 200, 0)]

# Scale the RGB values to the [0, 1] range
for i in range(len(colors)):
    r, g, b = colors[i]
    colors[i] = (r / 255., g / 255., b / 255.)

linew=2.0
with PdfPages(sys.argv[2]) as pdf:
	plt.figure(1)
	plt.title('Cells')
	plt.xlabel("Time")
	plt.ylabel("Quantity")
	plt.grid(True)
	for i, lab in enumerate(["uC", "RC", "WC", "DC", "WDC"]):
		plt.plot(data[:,0], data[:,i+1], c=colors[i], lw=linew, label=lab)
		if exp_data.get(lab):
			plt.plot(exp_times, exp_data[lab], "x", c=colors[i])
	if log_scale == 1:
		plt.yscale('log')
	plt.legend()
	pdf.savefig()
	plt.close()

	plt.figure(2)
	plt.title('Virions')
	plt.xlabel("Time")
	plt.ylabel("Quantity")
	plt.grid(True)
	for i, lab in enumerate(["WV", "DV"], 5):
		plt.plot(data[:,0], data[:,i+1], c=colors[i], lw=linew, label=lab)
		if exp_data.get(lab):
			plt.plot(exp_times, exp_data[lab], "x", c=colors[i])
	if log_scale == 1:
		plt.yscale('log')
	plt.legend()
	pdf.savefig()
	plt.close()


