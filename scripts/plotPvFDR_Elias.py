#! /usr/bin/env python
upFDR = 0.10
from pylab import  *
import glob

curves = []
curveNames = []
for doc in glob.glob('*.res'):
  print doc
  curveName = doc[:-4]
  truncName = False
  if curveName[-1]>="0" and curveName[-1]<="9":
    curveName = curveName[:-1]
    truncName=True
  fdrs = []
  passings = []
  fdr=.0
  tpATfdr=0
  passing,decoys,search,firstBreach=0,0,1,1
  f = open(doc,"r")
  labels = [int(line) for line in f.readlines()]
  f.close()
  for val in labels:
    passing += 1
    if val != 1:
      decoys+=1
    gfdr = 2*decoys/float(passing)
    if gfdr<upFDR:
      fdrs+=[gfdr]
      passings+=[passing]
  curves += [(truncName,(fdrs,passings),curveName)]
curves.sort()
curveNames = set()
for same,curve,curveName in curves:
  curveNames.add(curveName)
  if not same:
    plot(curve[0],curve[1])
  else:
    plot(curve[0],curve[1], 'g-')

figure(1).set_frameon(False)
xlabel('FDR',fontsize='large')
ylabel('Number of identifications',fontsize='large')
cn=list(curveNames)
#cn.reverse()
legend(cn,loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("fdrP_Elias.eps")
show()
