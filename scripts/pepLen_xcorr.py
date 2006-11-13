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
labels = [int(l.split()[1]) for l in ll[1:]]
xcorr=3
pepLen=7
pLen,nLen,px,nx=[],[],[],[]
for i in range(1,len(lf)):
  wf = [float(w) for w in lf[i].split()[1:]]
  if labels[i-1]==1:
    pLen += [wf[pepLen]]
    px += [wf[xcorr]]
  else:
    nLen += [wf[pepLen]]
    nx += [wf[xcorr]]
plot(pLen,px,'g+',label='normal')
plot(nLen,nx,'rx',label='shuffled')
xlabel('peptide length')
ylabel('Xcorr')
legend()
savefig('pepLen_Xcorr.eps')
savefig('pepLen_Xcorr.png')
show()
