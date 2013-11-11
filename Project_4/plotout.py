import sys
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from numpy import array, linspace
filenames=[]
files=[]
names=[]

if (len(sys.argv)<2):
	print "usage: python plotout.py filenames..." 
	sys.exit(1)
for name in sys.argv[1:]:
	filenames.append(name)
	files.append(open(name,'r'))
	names.append(name[-2:])
	

for i in range(len(files)):
		y=np.loadtxt(filenames[i])#,unpack=True)
		x=np.linspace(0,1,len(y[-1]))
		
		if((filenames[i])[-2:]=='EX'):
			plt.plot(x, y[-1]+1-x)
		else:
			plt.plot(x, y[-1]+1-x,'-')
		plt.xlabel('x')
		plt.ylabel('u')	
		plt.hold('on')
		#plt.ylim([-m,m])
		#plt.xlim([-m,m])
		plt.legend((names),prop={'size':8},labelspacing=0.25)
		
plt.savefig('%s.png' % 'infpntM')
plt.savefig('%s.pdf' % 'intpntM')
plt.figure()
for file in files:  file.close()

