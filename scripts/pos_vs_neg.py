#! /usr/bin/env python
import numpy
from pylab import  *
matplotlib.use('PS')
import glob
import re

x,y,z = [],[],[]
for doc in glob.glob('*.res'):
  print doc
  m = re.search('(-?[\d.]+)-(-?[\d.]+).res',doc)
  cpos=float(m.group(1))
  cneg=float(m.group(2))
#  curve = []
  fps = []
  tps = []
  fdr=.0
  tp,fp,search=0,0,1
  f = open(doc,"r")
  for line in f.readlines():
    val = int(line)
#    if search: print str(val) + " " + str(fdr)
    if val == 0:
      fp+=1
    if val == 1:
      tp+=1
      fdr = fp/(fp+float(tp))
#      curve += [fdr]
 #     print fdr
    fps+=[fp]
    tps+=[tp]
    if search == 1 and fdr > 0.01:
#      print curveName + ": FDR is 1% when finding " + str(tp) + " positives"
      break
  f.close()
  x += [cpos]
  y += [cneg]
  z += [tp]

xvec = list(set(x))
yvec = list(set(y))
xvec.sort()
yvec.sort()
X=numpy.array(xvec)
Y=numpy.array(yvec)
#X, Y = meshgrid(xvec,yvec)
Z=numpy.zeros((len(xvec),len(yvec)),'i')
#Z=zeros((len(xvec),len(yvec)))
for ix in range(len(z)):
  Z[xvec.index(x[ix])][yvec.index(y[ix])]=z[ix]
print Z
print X
ZZ=transpose(Z)
contourf(X,Y,ZZ)
xlabel('10log(Cpos)')
ylabel('10log(Cneg)')
#plot(fps,tps)
savefig("C.eps")
show()
