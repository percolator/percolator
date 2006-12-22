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
labels = [int(l.split()[1]) for l in ll[1:]]
xcorr=3
delt=2
tryC=11
tryN=12
pepSite=15
tp=0
fp=0
for i in range(1,len(lf)):
  wf = [float(w) for w in lf[i].split()[1:]]
  if(onlyCharge>0 and wf[onlyCharge]==0):
    continue
  if wf[delt]<0.08:
    continue
  if wf[pepSite]==0.0:
    continue
  if onlyTryp and (wf[tryC]!=1 or wf[tryN]!=1):
    continue
  if wf[8]==1:
    if wf[xcorr]<1.8:
      continue
  elif wf[9]==1:
    if wf[xcorr]<2.5:
      continue
  elif wf[10]==1:
    if wf[xcorr]<3.5:
      continue
  else:
    continue
  if labels[i]==1: tp+=1
  else: fp+=1
#  print wf[delt],wf[xcorr],wf[8],wf[9],wf[10],wf[tryN],wf[tryC],tp,fp
print tp, fp/(float(tp))
