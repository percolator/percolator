#! /usr/bin/env python
import sys
import math
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
peplen=7
pLen,nLen,px,nx=[],[],[],[]
for i in range(1,len(lf)):
  wf = [float(w) for w in lf[i].split()[1:]]
  pepLen = wf[peplen]
  Xcorr = wf[xcorr]
  if (wf[9]==1):
    # Charge 2
    if (pepLen>15):
      N=2*15
    else:
      N=2*pepLen
  elif (wf[10]==1):
    # Charge 3
    if (pepLen>25):
      N=4*25
    else:
      N=4*pepLen
  else:
    continue
  Xcorrp=math.log(Xcorr)/math.log(N)
  if labels[i-1]==1:
    pLen += [pepLen]
    px += [Xcorrp]
  else:
    nLen += [pepLen]
    nx += [Xcorrp]
plot(pLen,px,'g+',label='normal')
plot(nLen,nx,'rx',label='shuffled')
xlabel('peptide length')
ylabel('Xcorr prime')
legend()
savefig('pepLen_Xcorrp.eps')
savefig('pepLen_Xcorrp.png')
show()
