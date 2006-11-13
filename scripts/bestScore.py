#! /usr/bin/env python
import sys

def getScores(fileName,label):
  scores = []
  file = open(fileName,"r")
  bestScore=-1e100
  for line in file.readlines():
#  while (line=file.readline()):
    if (line[0]=='S'):
      if (bestScore>0): scores += [(bestScore,label)]
      bestScore=-1e100
    if (line[0]=='M'):
      score = float(line.split()[5])
      if (score>bestScore):
        bestScore=score
#      print score, bestScore
  return scores

scores=getScores(sys.argv[1],1)
scores+=getScores(sys.argv[2],-1)
scores.sort(reverse=True)
labels = [sc[1] for sc in scores]
for l in labels:
  print(str(l))
