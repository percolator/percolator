#! /usr/bin/env python
import sys
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
thresh = float(sys.argv[4])
onlyTryp = (len(sys.argv)>5 and sys.argv[5]=="Y")
onlyCharge = 0
if len(sys.argv)>6:
  onlyCharge = 7 + int(sys.argv[6])
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
names = [l.split()[0] for l in lf]
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
  scores += [(sum,label,names[i])]
scores.sort(reverse=True)
tp,fp,fdr=0,0,0
my_pi0=0.9
for sum,val,name in scores:
  if val != 1:
    fp+=1
  if val == 1:
    tp+=1
  if tp>0:
    fdr = my_pi0*fp/(float(tp))
  else:
    fdr = my_pi0
  if fdr<=thresh:
    if val==1: print name
  else:
    exit
