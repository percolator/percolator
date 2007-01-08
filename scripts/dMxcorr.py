#! /usr/bin/env python
import sys
from pylab import  *
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
lf=ffeatures.readlines()
ll=flab.readlines()
flab.close()
ffeatures.close()
scores = []
labels = [l.split()[1] for l in ll]
for i in range(1,len(lf)):
  sum = 0
  wf = lf[i].split()
  scores += [(wf[4],float(wf[2])/charge,int(labels[i]))]
plot([sc[1] for sc in scores if sc[2]==1],[sc[0] for sc in scores if sc[2]==1],'g+',label='Normal')
plot([sc[1] for sc in scores if sc[2]==-1],[sc[0] for sc in scores if sc[2]==-1],'rx',label='Shuffled')
xlabel("dM [Dalton]")
ylabel("Xcorr")
legend()
savefig("dMxcorrPlot.eps")
show()
