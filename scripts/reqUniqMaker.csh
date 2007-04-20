#!/usr/bin/env python
import sys

results=sys.stdin.readlines()
scoreDict=dict()
for line in results[1:]:
  wrd=line.split()
  k=wrd[3]
  q=float(wrd[2])
  if k in scoreDict:
    scoreDict[k]=min(scoreDict[k],q)
  else:
    scoreDict[k]=q
scores=scoreDict.values()
scores.sort()
for q in scores:
  print q