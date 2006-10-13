#! /usr/bin/env python
import sys

ffeatures = open(sys.argv[1],"r")
flab = open(sys.argv[2],"r")
fweights = open(sys.argv[3],"r")
wei=fweights.readline().split()
weights = [float(w) for w in wei]
lf=ffeatures.readlines()
ll=flab.readlines()
fweights.close()
flab.close()
ffeatures.close()
scores = []
labels = [l.split()[1] for l in ll]
for i in range(1,len(lf)):
  sum = 0
  wf = lf[i].split()
  for j in range(len(weights)-1):
    sum += float(wf[j+1])*weights[j]
  scores += [(sum,int(labels[i]))]
scores.sort(reverse=True)
labels = [sc[1] for sc in scores]
for l in labels:
  print(str(l))
