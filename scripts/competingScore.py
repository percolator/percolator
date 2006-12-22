#! /usr/bin/env python
import sys
ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
wl=fweights.readlines()
wei=wl[-1].split()
weights = [float(w) for w in wei]
lf=ffeatures.readlines()
ll=flab.readlines()
fweights.close()
flab.close()
ffeatures.close()
sdic = {}
scores = []
labels = [int(l.split()[1]) for l in ll[1:]]
for i in range(1,len(lf)):
  label=labels[i-1]
  if label<0: ixl= -label
  else: ixl=0
  sum = 0
  wf = lf[i].split()
  theId = wf[0]
  charno = theId.find('_')
  if charno >0:
    theId=theId[charno:]
  if not theId in sdic:
    sdic[theId]=(0,0,0)
  for j in range(len(weights)-1):
    sum += float(wf[j+1])*weights[j]
  sdic[theId][ixl]=sum
for theId in sdic.keys():
  if (sdic[theId][0]>sdic[theId][1]):
    scores += [(sdic[theId][0],1)]
  if (sdic[theId][2]>sdic[theId][1]):
    scores += [(sdic[theId][0],-1)]
scores = sdic.values()
scores.sort(reverse=True)
labels = [sc[1] for sc in scores]
for l in labels:
  print(str(l))
