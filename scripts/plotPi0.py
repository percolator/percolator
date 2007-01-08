#! /usr/bin/env python
upFDR = 0.10
from pylab import  *
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

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
    if lam<0.95:
      fdrs+=[lam]
      pi0+=[min(1,(1-plam)/(1-lam))]
  curves += [((fdrs, pi0),curveName)]
for curve,curveName in curves:
  plot(curve[0],curve[1],label=curveName)
xlabel('Lambda',fontsize='large')
ylabel('Pi0',fontsize='large')
legend(loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("pi0.eps")
show()
