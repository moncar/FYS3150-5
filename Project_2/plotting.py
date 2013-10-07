from scitools.std import *
import numpy as np

#Reading vector:
def vecread(filename):
	table = open(filename, 'r')
	rows = table.readlines()
	n = 1
	m = len(rows)
	#print m
	c = 0
	elements = np.zeros(m)
	for i in rows:
		row = i.split(' ')
		row = ''.join(row)
		#row = filter(lambda number: number.strip(), row)
		#print row
		for k in range(len(row)):
			elements[c] = float(row)
		c += 1
	return elements

u1 = vecread("output_0.01")
u2 = vecread("output_0.5")
u3 = vecread("output_1")
u4 = vecread("output_5")
#print len(v), len(vluf)
x = linspace(0,5,len(u1))
#u = 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x)

#Making and saving plots:
figure(1)
plot(x,u1, xlabel = "rho", ylabel = "u(x)", legend = "omega_r = 0.01")
hold('on')
plot(x,u2**2.0, legend = "omega_r = 0.5")
plot(x,u3**2.0, legend = "omega_r = 1")
plot(x,u4**2.0, legend = "omega_r = 5",hardcopy = "waveplot.eps")
