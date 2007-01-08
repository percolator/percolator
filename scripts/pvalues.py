#! /usr/bin/env python
import glob

def countNeg(x,y):
  if (y!=1):
    x+=1
  return x

curves = []
curveNames = []
for doc in glob.glob('*.res'):
  print doc
  outName = doc[:-4] + ".pval"
  f = open(doc,"r")
  labels = [int(line) for line in f.readlines()]
  f.close()
  fp,tp=0,0
  f = open(outName,"w")
  nneg = max(float(reduce(countNeg,labels,0)),69705.0)
  npos = max(float(len(labels) - nneg),69705.0)
  for val in labels:
    if val != 1:
      fp+=1
    pval=fp/nneg
    if val == 1:
      f.write(str(pval))
      f.write('\n')
  f.close()
