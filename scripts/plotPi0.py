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
  nneg = float(reduce(countNeg,labels,0))
  npos = float(len(labels) - nneg)
  found = False
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
    if lam<0.99:
      fdrs+=[lam]
      pi0+=[min(1,(1-plam)/(1-lam))]
    if (not found and lam > 0.10):
      found = True
      level = 1 - plam
  curves += [(level,(fdrs, pi0),curveName)]
curves.sort(reverse=True)
for level,curve,curveName in curves:
  plot(curve[0],curve[1],label=curveName)
xlabel('Lambda',fontsize='large')
ylabel('Pi0',fontsize='large')
legend(loc=3)
savefig("pi0.png")
savefig("pi0.eps")
show()
