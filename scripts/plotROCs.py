#! /usr/bin/env python
upFDR = 0.10
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
    if val != 1:
      fp+=1
    if val == 1:
      tp+=1
      fdr = fp/(float(tp))
    if fdr>upFDR:
      continue
    fdrs+=[fdr]
    tps+=[tp]
    if search == 1 and fdr > 0.01:
      print curveName + ": FDR is 1% when finding " + str(tp) + " positives"
      search = 0
      tpATfdr = tp
  f.close()
  curves += [(tpATfdr,(fdrs, tps),curveName)]
curves.sort(reverse=True)
for tp,curve,curveName in curves:
  plot(curve[0],curve[1],label=curveName)
style = "o^sDvx+<>"
i = 0
for doc in glob.glob('*.pnt'):
  print doc
  curveName = doc[:-4]
  f = open(doc,"r")
  line = f.readline()
  f.close
  field = [float(w) for w in line.split()]
  plot([field[1]],[field[0]],style[i],label=curveName)
  i += 1

#axis([0,0.05,tp-200,tp+300])
#axis([0,0.05,0,16000])
#axis([0,0.05,0,200])
xlabel('False Discovery Rate')
ylabel('True Positives')
legend(loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02)
savefig("roc.eps")
show()
