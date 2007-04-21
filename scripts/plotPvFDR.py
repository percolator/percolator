#! /usr/bin/env python
upFDR = 0.10
from pylab import  *
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

curves = []
curveNames = set()
pi0 = 1.0
for doc in glob.glob('*.res'):
  print doc
  curveName = doc[:-4]
  truncName = False
  if curveName[-1]>="0" and curveName[-1]<="9":
    curveName = curveName[:-1]
    truncName=True
  fdrs = []
  tps = []
  fdr=.0
  tpATfdr=0
  tp,fp,search,firstBreach=0,0,1,1
  f = open(doc,"r")
  labels = [int(line) for line in f.readlines()]
  f.close()
  if curveName.find("rank")>0:
    my_pi0 = 1 - (1-pi0)/5.0
    print "pi0 div with 5"
  else:
    my_pi0=pi0
  for val in labels:
    if val != 1:
      fp+=1
    if val == 1:
      tp+=1
    if tp>0:
      fdr = my_pi0*fp/(float(tp))
    else:
      fdr = 1
#    if fdr>upFDR:
#      continue
    if fdr<upFDR:
      fdrs+=[fdr]
      tps+=[tp]
  for i in range(len(fdrs)-1,0,-1):
    if search == 1 and fdrs[i-1] <= 0.01:
      print curveName + ": FDR is 1% when finding " + str(tps[i-1]) + " positives"
      search = 0
      tpATfdr = tps[i-1]
  curves += [(truncName,tpATfdr,(fdrs, tps))]
  curveNames.add(curveName)

curveNames=list(curveNames)
#curves.sort(reverse=True)
curves.sort()
for same,tp,curve in curves:
  if not same:
    plot(curve[0],curve[1])
  else:
    plot(curve[0],curve[1], 'g-')

figure(1).set_frameon(False)
xlabel('FDR',fontsize='large')
ylabel('Number of identifications',fontsize='large')
legend(curveNames,loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("pVsFDR.eps")
show()
