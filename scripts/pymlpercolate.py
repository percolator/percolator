#! /usr/bin/env python
import numpy
import datafunc
import svm

def scoreLinear(X,lab,w,fn):
  scores=[]
  for ix in range(len(X)):
    score = 0.0
    for i in range(len(w)):
      score += X[ix][i]*w[i]
    scores += [(score,lab[ix])]
  scores.sort(reverse=True)
  labels = [sc[1] for sc in scores]
  f=open(fn + '.res',"w")
  for l in labels:
    f.write(str(l)+'\n')
  f.close

def getTrainingData(X,lab,w):
  scores=[]
  keep_ixs=[]
  for ix in range(len(X)):
    score = 0.0
    for i in range(len(w)):
      score += X[ix][i]*w[i]
      scores += [(score,ix,lab[ix])]
      if lab[ix]==0:
        keep_ixs += [ix]
  scores.sort(reverse=True)
  fp,tp,fdr=0,0,0
  for score,ix,l in scores:
    if l == 0:
      fp+=1
    else:
      tp+=1
    fdr = fp/(fp+float(tp))
    if fdr>0.01:
      break
    if l == 1:
      keep_ixs += [ix]
  XX = numpy.array([X[ix] for ix in keep_ixs])
  ll = numpy.array([lab[ix] for ix in keep_ixs])
  return datafunc.DataSet(X=XX,L=ll)

wholedata = datafunc.DataSet("6d.pyml", labelsColumn = -1)
w=numpy.array([0,1,0,0])
for it in range(3):
  d=getTrainingData(wholedata.X,wholedata.labels.Y,w)
  s = svm.SVM()
#  s.Cmode="equal"
  s.train(d)
  w=s.model.w
print("Model parameters:")
print(s.model.b)
print(s.model.w)
fn = "Itt-fdr-6d"
scoreLinear(wholedata.X,wholedata.labels.Y,s.model.w,fn)
