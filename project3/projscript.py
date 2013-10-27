import os, sys

usage="""
python  projscript.py x [ID] [filename] [Nk]
x =0 : no run
x =1 : run
x =2 : plot to file 
ID=0 : Sun Earth
ID=2 : All planets + Sun

"""
if len(sys.argv)<2:
	print usage
	sys.exit(1)

execute=int(sys.argv[1])
ID='0 '
fn='fileout'
Nk=' 1 '	
comp='g++ '
files = ['project3','celeb','planetarysystem']
py='python '


if len(sys.argv)==5:
	ID=sys.argv[2]+' '
	fn=sys.argv[3]
	Nk=' '+sys.argv[4]+' '


ocmd=comp+'-o run.x '
for f in files:#compile
	cmd = comp+'-c '+f+'.cpp' 
	print cmd
	failure = os.system(cmd) 
	if failure:
		print "Command failed", cmd;sys.exit(1)
	ocmd +=f+'.o '
cmd=ocmd
print cmd #link
failure = os.system(cmd) 
if failure:
	print "Command failed", cmd;sys.exit(1)
if execute!=0:#execute
	cmd='./run.x '+ID+fn+Nk
	print cmd
	failure = os.system(cmd) 
	if failure:
		print "Command failed", cmd;sys.exit(1)
if execute>=2:#plot
	cmd=py+'plotout.py '+ID+' '+fn 
	print cmd
	failure = os.system(cmd) 
	if failure:
		print "Command failed", cmd;sys.exit(1)
if execute==3:#show
	cmd = 'xdg-open '+fn+'.pdf'
	print cmd
	failure = os.system(cmd)
if execute==-1:#readfile
	cmd = 'cat '+fn
	print cmd
	failure = os.system(cmd)

