#! /usr/bin/env python
from pylab import  *
matplotlib.use('PS')
import glob

curves = []
curveNames = []
for doc in glob.glob('*.res'):
  print doc
  curveName = doc[:-4]
#  curve = []
#  fps = []
  fdrs = []
  tps = []
  fdr=.0
  tpATfdr=0
  tp,fp,search=0,0,1
  f = open(doc,"r")
  for line in f.readlines():
    val = int(line)
#    if search: print str(val) + " " + str(fdr)
    if val != 1:
      fp+=1
    if val == 1:
      tp+=1
      fdr = fp/(fp+float(tp))
#      curve += [fdr]
 #     print fdr
#    fps+=[fp]
    fdrs+=[fdr]
    tps+=[tp]
    if search == 1 and fdr > 0.01:
      print curveName + ": FDR is 1% when finding " + str(tp) + " positives"
      search = 0
      tpATfdr = tp
  f.close()
#  curves += [([fp/float(fps[-1]) for fp in fps],[tp/float(tps[-1]) for tp in tps])]
  curves += [(tpATfdr,(fdrs, tps),curveName)]
#  curveNames += [curveName]
curves.sort(reverse=True)
#a=axes([0,0.1,0,0.15])
#for ix in range(len(curveNames)):
for tp,curve,curveName in curves:
  plot(curve[0],curve[1],label=curveName)
axis([0,0.03,400,800])
xlabel('False Discovery Rate')
ylabel('True Positives')
legend()
#plot(fps,tps)
savefig("roc.eps")
show()
