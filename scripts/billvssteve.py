#! /usr/bin/env python
from pylab import  *

def readFdrs(fn):
  fil=open(fn,"r")
  lines=fil.readlines()
  fil.close()
  valarr = [(float(v[0]),float(v[1])) for v in [l.split() for l in lines]]
  return valarr

def findFirst(tplev,arr):
  for tp,fdr in arr:
    if tplev<tp: return fdr
  return 1.0

steve = readFdrs(sys.argv[1])
bill = readFdrs(sys.argv[2])

sfdrs=[]
bfdrs=[]
for tp,bfdr in bill:
  sfdr = findFirst(tp,steve)
  bfdrs += [bfdr]
  sfdrs += [sfdr]

plot(bfdrs,sfdrs)
plot([0,bfdrs[-1]],[0,bfdrs[-1]])

xlabel('Bill FDR')
ylabel('Steve FDR')
savefig("billvssteve.eps")
show()
