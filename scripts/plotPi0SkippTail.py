#!/usr/bin/env python
upFDR = 0.10
from pylab import  *
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

tail=0.4
curves = []
curveNames = []
for doc in glob.glob('*.res'):
  print doc
  curveName = doc[:-4]
#  curve = []
#  fps = []
  fdrs = []
  pi0 = []
  fdr,oldfdr=.0,.0
  tpATfdr=0
  tp,fp,search,firstBreach=0,0,1,1
  f = open(doc,"r")
  labels = [int(line) for line in f.readlines()]
  f.close()
  nneg = max(float(reduce(countNeg,labels,0)),69705.0)
  npos = max(float(len(labels) - nneg),69705.0)
  plamt=-1
  for val in labels:
    if val != 1:
      fp+=1
    if val == 1:
      tp+=1
    lam = fp/nneg
    plam = tp/npos
    if lam>=tail and plamt<0:
      plamt=plam
  print plamt
  tp,fp=0,0
  for val in labels:
    if val != 1:
      fp+=1
    if val == 1:
      tp+=1
    lam = fp/nneg
    plam = tp/npos
    if lam<tail-0.1:
      fdrs+=[lam]
      pi0+=[min(1,(plamt-plam)/(tail-lam)*(1-tail)/(1-plamt))]
  if plamt>0:
    curves += [((fdrs, pi0),curveName)]
for curve,curveName in curves:
  plot(curve[0],curve[1],label=curveName)
xlabel('Lambda',fontsize='large')
ylabel('Pi0',fontsize='large')
legend(loc=1,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("pi0.png")
savefig("pi0.eps")
show()
