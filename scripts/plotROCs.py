#! /usr/bin/env python
upFDR = 0.10 - 1e-5
orderFile = "order"
from pylab import  *
#from matplotlib import rc
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
print rc
curves = []
curveNames = []
pi0 = 0.9
if glob.glob(orderFile):
  f = open(orderFile,"r")
  curveFiles = [line[:-1] for line in f.readlines()]
  f.close()
  preSorted=True
else:
  curveFiles = glob.glob('*.re[qs]')
  preSorted=False

for doc in curveFiles:
  print doc
  curveName = doc[:-4]
#  curve = []
#  decoys = []
  fdrs = []
  targets = []
  fdr=.0
  targetATfdr=0
  target,decoy,search,firstBreach=0,0,1,1
  if doc[-1]=="q":
    # req extension
    f = open(doc,"r")
    fdrs = [float(line) for line in f.readlines()]
    f.close()
    fdrs = [fdr for fdr in fdrs if fdr<upFDR]
    fdrs.sort()
    targets=range(1,len(fdrs)+1)
  else:
    # res extension
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
        decoy+=1
      if val == 1:
        target+=1
      if target>0:
        fdr = my_pi0*decoy/(float(target))
      else:
        fdr = 1
      if fdr<upFDR:
        fdrs+=[fdr]
        targets+=[target]
    if len(targets)==0:
      fdrs+=[0]
      targets+=[0]
    fdrs+=[upFDR]
    targets+=[targets[-1]]
  for i in range(len(fdrs)-1,0,-1): # Make q-values out of the fdrs
    if fdrs[i-1]>fdrs[i]: fdrs[i-1]=fdrs[i]
    if search == 1 and fdrs[i-1] <= 0.01:
      print curveName + ": FDR is 1% when finding " + str(targets[i-1]) + " positives"
      search = 0
      targetATfdr = targets[i-1]
  curves += [(targetATfdr,(fdrs, targets),curveName)]

if not preSorted:
  curves.sort(reverse=True)
for target,curve,curveName in curves:
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
xlabel('q value',fontsize='large')
ylabel('Number of peptide-spectrum matches identified',fontsize='large')
legend(loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02,numpoints=3)
savefig("roc.eps")
show()
