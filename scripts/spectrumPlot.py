#! /usr/bin/env python
import sys
from pylab import  *

def plotSpec(linesPeaks,key,color):
  colPeaks=[(line.split()[0],line.split()[1]) for line in linesPeaks if line.count(key)]
  if len(colPeaks)==0: return
  colSpectra=[]
  for peak in colPeaks:
    colSpectra += [(peak[0],0)]
    colSpectra += [peak]
    colSpectra += [(peak[0],0)]
  plot([peak[0] for peak in colSpectra],[peak[1] for peak in colSpectra],color)

fhPeaks = open(sys.argv[1],"r")
linesPeaks=fhPeaks.readlines()
if linesPeaks[0][0]=='>':
  linesPeaks=linesPeaks[1:]
fhPeaks.close()
plotSpec(linesPeaks,"",'y')
plotSpec(linesPeaks,"magenta",'m')
plotSpec(linesPeaks,"red",'r')
plotSpec(linesPeaks,"blue",'b')
xlabel('m/z')
ylabel('Abundance')
#legend(loc=4,pad=0.1,labelsep = 0.001,handlelen=0.04,handletextsep=0.02)
savefig("spectrum.eps")
show()
