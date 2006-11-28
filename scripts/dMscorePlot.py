#! /usr/bin/env python
import sys
from pylab import  *
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
wl=fweights.readlines()
wei=wl[-1].split()
weights = [float(w) for w in wei]
lf=ffeatures.readlines()
ll=flab.readlines()
fweights.close()
flab.close()
ffeatures.close()
scores = []
labels = [l.split()[1] for l in ll]
for i in range(1,len(lf)):
  sum = 0
  wf = lf[i].split()
  for j in range(len(weights)-1):
    sum += float(wf[j+1])*weights[j]
  if (sum>5):
    scores += [(sum,wf[2],int(labels[i]))]
plot([sc[1] for sc in scores if sc[2]==1],[sc[0] for sc in scores if sc[2]==1],'g+',label='Normal')
plot([sc[1] for sc in scores if sc[2]==-1],[sc[0] for sc in scores if sc[2]==-1],'rx',label='Shuffled')
savefig("dMscorePlot.eps")
xlabel("dM [Dalton]")
ylabel("Percolator score")
legend()
show()