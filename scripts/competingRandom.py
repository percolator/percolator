#! /usr/bin/env python
# script that emaulates the effeects of a combined search.
# the script sorts scores based on the bigger of the scores
# from normal and testset, sorts the scores an report a 1 if
# the normal hit was highest -1 otherwhise.
#
downSample=100
tv=5.0
fv=1.5
pi0=0.9
import sys
import os
import random

ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
#onlyTryp = (len(sys.argv)>4 and sys.argv[4]=="Y")
wl=fweights.readlines()
wei=wl[-1].split()
weights = [float(w) for w in wei]
lf=ffeatures.readlines()
ll=flab.readlines()
fweights.close()
flab.close()
ffeatures.close()
sdic = {}
tryp = {}
scores = []
labels = [int(l.split()[1]) for l in ll[1:]]

for i in range(1,len(lf)):
  label=labels[i-1]
  if label<0: ixl= -label
  else: ixl=0
  sum = 0
  wf = lf[i].split()
  if float(wf[10]) != 1.0:
    continue
  theId = wf[0]
  charno = theId.find('_')
  if charno >0:
    theId=theId[charno:]
  if not theId in sdic:
    sdic[theId]=[0,0,0]
    tryp[theId]=[False,False,False]
  for j in range(len(weights)-1):
    sum += float(wf[j+1])*weights[j]
  sdic[theId][ixl]=sum
  tryp[theId][ixl]=(float(wf[12])>0) and (float(wf[13])>0)
decoyScores=[v[2] for k,v in sdic.items()]
#decoyScores.sort(reverse=True)
random.shuffle(decoyScores)
specN,randN,ix=0,0,0
spectra=float(len(sdic))
for theId in sdic.keys():
  if sdic[theId][0]>sdic[theId][2]:
    specN += 1
  if sdic[theId][0]>decoyScores[ix]:
    randN += 1
  ix+=1
print specN,specN/spectra,randN,randN/spectra,spectra
