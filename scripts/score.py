#! /usr/bin/env python
import sys
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
onlyTryp = (len(sys.argv)>4 and sys.argv[4]=="Y")
onlyCharge = 0
if len(sys.argv)>5:
  onlyCharge = 7 + int(sys.argv[5])
wl=fweights.readlines()
wei=wl[-1].split()
weights = [float(w) for w in wei]
lf=ffeatures.readlines()
ll=flab.readlines()
fweights.close()
flab.close()
ffeatures.close()
scores = []
labels = [l.split()[1] for l in ll]
for i in range(1,len(lf)):
  label = int(labels[i])
  if (label==-1): continue
  if (label==-2): label=-1
  sum = 0
  wf = [float(w) for w in lf[i].split()[1:]]
  if(onlyCharge>0 and wf[onlyCharge]==0):
    continue
  if (onlyTryp and ((wf[11]==0) or (wf[12]==0))):
    continue
  for j in range(len(weights)-1):
    sum += float(wf[j])*weights[j]
  scores += [(sum,label)]
scores.sort(reverse=True)
labels = [sc[1] for sc in scores]
for l in labels:
  print(str(l))
