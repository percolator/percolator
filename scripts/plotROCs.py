#! /usr/bin/env python
upFDR = 0.10 - 1e-5
orderFile = "order"
from pylab import  *
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

curves = []
curveNames = []
pi0 = 0.9
if glob.glob(orderFile):
  f = open(orderFile,"r")
  curveFiles = [line[:-1] for line in f.readlines()]
  f.close()
  preSorted=True
else:
  curveFiles = glob.glob('*.res')
  preSorted=False

for doc in curveFiles:
  print doc
  curveName = doc[:-4]
#  curve = []
#  fps = []
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
  if len(tps)==0:
    fdrs+=[0]
    tps+=[0]
  fdrs+=[upFDR]
  tps+=[tps[-1]]
  for i in range(len(fdrs)-1,0,-1): # Make q-values out of the fdrs
    if fdrs[i-1]>fdrs[i]: fdrs[i-1]=fdrs[i]
    if search == 1 and fdrs[i-1] <= 0.01:
      print curveName + ": FDR is 1% when finding " + str(tps[i-1]) + " positives"
      search = 0
      tpATfdr = tps[i-1]
  curves += [(tpATfdr,(fdrs, tps),curveName)]

if not preSorted:
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
  if field[1]>upFDR:
    continue
  plot([field[1]*pi0],[field[0]],style[i],label=curveName)
  i += 1

figure(1).set_frameon(False)
xlabel('q-value',fontsize='large')
ylabel('Number of peptide-spectrum matches identified',fontsize='large')
legend(loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("roc.eps")
show()
