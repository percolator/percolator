#! /usr/bin/env python
import numpy
from pylab import  *
matplotlib.use('PS')
import math
import glob
import re

x,y,z = [],[],[]
hightp,hx,hy = {},{},{}
for doc in glob.glob('*.res'):
#  print doc
  m = re.search('([\d.]+)-([\d.]+)-([\d.]+).res',doc)
  FDR=float(m.group(1))
  cpos=float(m.group(2))
  cneg=float(m.group(3))
  fdr=.0
  tp,fp,search=0,0,1
  f = open(doc,"r")
  for line in f.readlines():
    val = int(line)
#    if search: print str(val) + " " + str(fdr)
    if val == 1:
      tp+=1
    else:
      fp+=1
    fdr = fp/(fp+float(tp))
#      curve += [fdr]
 #     print fdr
    if fdr > 0.01:
#      print doc + ": FDR is 1% when finding " + str(tp) + " positives"
      break
#    if (doc=="0.01-100-100.res"):
#      print line[:-1] + " " + str(val) + " " + str(tp) + " " + str(fp) + " " + str(fdr)
  f.close()
  if FDR == 0.03:
#  if cpos == cneg:
    x += [math.log(cpos,10)]
    y += [math.log(cneg,10)]
    z += [tp]
  ff = str(FDR)
  if (not hightp.has_key(ff)) or hightp[ff]<tp:
    hightp[ff]=tp
    hx[ff]=cpos
    hy[ff]=cneg

print hightp
print hx
print hy

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
#print Z
#print X
ZZ=transpose(Z)
CS=contourf(X,Y,ZZ,255)
xlabel('10log(Cpos)')
ylabel('10log(Cneg)')
jet()
colorbar(clabels=('hi','lo'))
#ylabel('FDR')
#plot(fps,tps)
savefig("C.eps")
show()
