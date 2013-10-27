import sys
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from numpy import array, linspace
filenames=[]
files=[]
names=[]


#usage: python plotout.py case filenames... 
c=["Sun","Earth","Moon","Mercury","Venus","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto",'COM']
case = range(5)
for i in case:
	if int(sys.argv[1])==i:
		case = i
resu = [3,4,11,4,3]
celebs=[
	[c[0],c[1],c[-1]],
	[c[0],c[1],c[2],c[-1]],
	[c[0],c[1],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[-1]],
	[c[0],c[1],c[7],c[-1]],
	[c[0],c[3],c[-1]]
       ]
mval=[2, 2, 40, 20, 0.5]

for name in sys.argv[2:]:
	filenames.append(name)
	files.append(open(name,'r'))
	names.append(name.split('.')[0])

length=2*resu[case]
m=mval[case]
for i in range(len(files)):
	col=0
	while (col<length):
		cols=np.loadtxt(filenames[i], usecols=(col,col+1),unpack=True)
		#plot col 1 vs col 2
		if col==length-2:
			plt.plot(cols[0], cols[1],'o',markersize=3)
		elif col==0:
			plt.plot(cols[0], cols[1],'o')
		else:
			plt.plot(cols[0], cols[1],'-')
		plt.xlabel('x[AU]')
		plt.ylabel('y[AU]')	
		plt.hold('on')
		plt.ylim([-m,m])
		plt.xlim([-m,m])
		#plot col 1 & col 2 vs x
		"""
		x=linspace(0,1,shape(cols)[0])		
		plt.plot(x, v1,label='v(x)')
		plt.hold('on')
		plt.plot(x, v2,'--o',label='u(x)')
		"""
		col+=2
	plt.legend(celebs[case],prop={'size':8},labelspacing=0.25)
		
	plt.savefig('%s.png' % names[i])
	plt.savefig('%s.pdf' % names[i])
	plt.figure()

for file in files:  file.close()

