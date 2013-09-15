import sys
import matplotlib.pyplot as plt
from math import log10
from numpy import array, linspace
filenames=[]
files=[]

for name in sys.argv[1:]:
	filenames.append(name)
	files.append(open(name,'r'))

for i in range(len(files)):
	values=[]
	for line in files[i].readlines():
		values.append([float(n) for n in line.strip().split(' ')])
		
	v1,v2 = zip(*values)
	v1=array(v1) ; v2=array(v2)
		
	x=linspace(0,1,len(values))
	plt.plot(x, v1,label='v(x)')
	plt.hold('on')
	plt.plot(x, v2,'--o',label='u(x)')
	plt.legend()
	plt.xlabel('x')
	plt.ylabel('f(x)')
		
	plt.title('log(n)=%.3g' %(log10(len(values)-2)))
	plt.savefig('poisson_%g.png' % (int(log10(len(values)))))
	plt.savefig('poisson_%g.pdf' % (int(log10(len(values)))))
	plt.figure()
for file in files:  file.close()
