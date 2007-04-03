#! /usr/bin/env python
import sys
import math
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
onlyTryp = (len(sys.argv)>3 and sys.argv[3]=="Y")
onlyCharge = 0
if len(sys.argv)>4:
  onlyCharge = 7 + int(sys.argv[4])
lf=ffeatures.readlines()
ll=flab.readlines()
flab.close()
ffeatures.close()
scores = []
labels = [l.split()[1] for l in ll]
for i in range(1,len(lf)):
  sum = 0
  wf = [float(w) for w in lf[i].split()[1:]]
  label = int(labels[i])
  if (label==-1): # only do test labels
    continue
  if (label==-2):
    label=-1
  if(onlyCharge>0 and wf[onlyCharge]==0):
    continue
  if (onlyTryp and ((wf[11]==0) or (wf[12]==0))):
    continue
  rSp = math.log(float(wf[0]))
  dM = abs(float(wf[1]))
  dCn = float(wf[2])
  Xcorr = float(wf[3])
  pepLen = int(wf[7])
  if (Xcorr<=0.0):
    Xcorr = 1e-10
  if (wf[9]==1):
    # Charge 2
    if (pepLen>15):
      N=2*15
    else:
      N=2*pepLen
    Xcorrp=math.log(Xcorr)/math.log(N)
    score = 8.362*Xcorrp + \
            7.386*dCn + \
            -0.194*rSp + \
            -0.314*dM + \
            -0.959
    scores += [(score,label)]
  elif (wf[10]==1):
    # Charge 3
    if (pepLen>25):
      N=4*25
    else:
      N=4*pepLen
#    print str(Xcorr) + " " + str(N)
    Xcorrp=math.log(Xcorr)/math.log(N)
    score = 9.933*Xcorrp + \
            11.149*dCn + \
            -0.201*rSp + \
            -0.277*dM + \
            -0.1460
    scores += [(score,labels[i])]
scores.sort(reverse=True)
labels = [sc[1] for sc in scores]
for l in labels:
  print(str(l))
