#! /usr/bin/env python

import Numeric
from pylab import  *
import math

matplotlib.use('PS')

def getWeights(fileName):
  fweights = open(fileName,"r")
  wl=fweights.readlines()
  fweights.close()
  wei=wl[-1].split()
  weights = [float(w) for w in wei]
  return weights

def getThreshold(scores):
  scores.sort(reverse=True)
  targets,decoys=0.0,0.0
  for sc,label in scores:
    if label>0:
      targets+=1
    else:
      decoys+=1
    if decoys/targets*0.9>0.01:
      return sc

w1 = getWeights(sys.argv[1])
w2 = getWeights(sys.argv[2])

ffeatures = open(sys.argv[3],"r")
flab = open(sys.argv[4],"r")
lf=ffeatures.readlines()
ll=flab.readlines()
flab.close()
ffeatures.close()

xtext = sys.argv[5] + " scoring function"
ytext = sys.argv[6] + " scoring function"

cloud = []
xx,yy,xxscores,yyscores = [],[],[],[],[]
labels = [l.split()[1] for l in ll]

for i in range(1,len(lf)):
  label = int(labels[i])
  if label==-1: continue
  sum1 = w1[-1]
  sum2 = w2[-1]
  wf = [float(w) for w in lf[i].split()[1:]]
  for j in range(len(wf)):
    sum1 += wf[j]*w1[j]
    sum2 += wf[j]*w2[j]
  cloud +=[(sum2,sum1)]
  xx +=[sum1]
  yy +=[sum2]
  xxscores += [(sum1,label)]
  yyscores += [(sum2,label)]

xt=getThreshold(xxscores)
yt=getThreshold(yyscores)

nX, mybinsX, patches = hist(xx, bins = 100)
#plot(mybins)
savefig("xx.eps")
clf()
nd, mybinsd, patches = hist(yy, bins = 100)
savefig("yy.eps")
clf()
cnt = Numeric.zeros((len(mybinsd),len(mybinsX)),'i')
#lcnt = Numeric.zeros((len(mybinsd),len(mybinsX)),'f')
for pt in cloud:
  d,x = 0,0
  while mybinsX[x] <= pt[1] and x + 1 < len(mybinsX):
    x+=1
  while mybinsd[d] <= pt[0] and d + 1 < len(mybinsd):
    d+=1
  cnt[d][x] += 1
a = [x for x in mybinsX]
#a += [a[-1]+.01,7]
mybinsX = array(a)
a = [x for x in mybinsd]
#a += [a[-1]+.01,0.7]
mybinsd = array(a)
#for a in range(len(mybinsd)):
#  for b in range(len(mybinsX)):
#    lcnt[a][b] = log(.000001+cnt[a][b])
#print lcnt
jet()
contourf(mybinsX,mybinsd,cnt,20)
#colorbar()
axvline(x=xt,color='k')
axhline(y=yt,color='k')
xlabel(xtext,fontsize='large')
ylabel(ytext,fontsize='large')
savefig("cloud.eps")
clf()
show()
#
# left, down = 0.3, 0.7
# rect1 = [0.1, 0.3, 0.2, 0.6]
# rect2 = [0.3, 0.3, 0.6, 0.6]
# rect3 = [0.3, 0.1, 0.6, 0.2]
#
# axBottom      = axes(rect3)#,sharex=axMiddle)
# a = [x for x in nX]
# a += [0,0]
# nX=array(a)
# axBottom.bar(mybinsX,nX,width=.01)
# 
# axMiddle = axes(rect2,sharex=axBottom)
# axMiddle.contourf(mybinsX,mybinsd,cnt,20)
# 
# a = [x for x in nd]
# a += [0,0]
# nd=array(a)
# 
# axLeft       = axes(rect1,sharey=axMiddle)  #left, bottom, width, height
# axLeft.barh(nd,mybinsd,height=.0)
# 
# #axBottom.xlim(axMiddle.xlim())
# 
# axLeft.set_ylabel('Trypsin')
# axBottom.set_xlabel('Elastase')
# setp( axLeft.get_xticklabels(), visible=False)
# setp( axBottom.get_yticklabels(), visible=False)
# setp( axMiddle.get_xticklabels(), visible=False)
# setp( axMiddle.get_yticklabels(), visible=False)
# #setp(axLeft, xticks=[])
# #setp(axBottom, yticks=[])
# 
# savefig("all.eps")
# show()
# 
# 
