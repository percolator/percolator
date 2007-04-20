#!/usr/bin/env python
import sys

results=sys.stdin.readlines()
scoreDict=dict()
scores = []
for line in results[1:]:
  wrd=line.split()
  scores+=[float(wrd[2])]
scores.sort()
for q in scores:
  print q
