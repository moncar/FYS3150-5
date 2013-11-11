#generic compilation of files(+run) with args
import os, sys
usage='usage: python cmparg.py execute filename args..'
if (len(sys.argv)<3):
	print usage
	sys.exit(0)
f=sys.argv[2]
args=sys.argv[3:]
x=int(sys.argv[1])
comp='g++ '
cmd = comp+'-c '+f+'.cpp' 
ocmd=comp+'-o '+f+'.x '+f+'.o'

print cmd
failure = os.system(cmd) 
if failure:
	print "Command failed", cmd;sys.exit(1)

cmd=ocmd
print cmd #link
failure = os.system(cmd) 
if failure:
	print "Command failed", cmd;sys.exit(1)

if x!=0:
	cmd='./'+f+'.x'
	for arg in args:
		cmd+=' '+arg
	print cmd
	failure = os.system(cmd) 
	if failure:
		print "Command failed", cmd;sys.exit(1)
