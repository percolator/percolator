#! /usr/bin/env python
upFDR = 0.10
from pylab import  *
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

pv = []
rankpv = []
fdr,oldfdr=.0,.0
tpATfdr=0
tp,fp,search,firstBreach=0,0,1,1
f = open(sys.argv[1],"r")
labels = [int(line) for line in f.readlines()]
f.close()
nneg = max(float(reduce(countNeg,labels,0)),69705.0)
npos = max(float(len(labels) - nneg),69705.0)
for val in labels:
  if val != 1:
    fp+=1
  if val == 1:
    tp+=1
  if tp>0:
    fdr = fp/(float(tp))
  else:
    fdr = 1
  if fdr<oldfdr:
    fdr=oldfdr
  oldfdr=fdr
#    if fdr>upFDR:
#      continue
  lam = fp/nneg
  plam = tp/npos
  pv +=[lam]
  rankpv+=[plam]
loglog([0.0001,1],[0.0001,1])
loglog(pv,rankpv)
#plot([0.0001,1],[0.0001,1])
#plot(pv,rankpv)
xlabel('empirical p-value',fontsize='large')
ylabel('Rank of empirical p-value',fontsize='large')
#savefig("pi0.eps")
show()
