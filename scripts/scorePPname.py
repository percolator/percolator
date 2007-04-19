#!/usr/bin/env python
import sys
fprob = open(sys.argv[1],"r")
lprob=fprob.readlines()
fprob.close()

fesi = open(sys.argv[2],"r")
lesi=fesi.readlines()
fesi.close()

s=[(line.split()[0],(float(line.split()[1]),line.split()[0].find("shuffled")<0,"","")) for line in lprob]
scoreDict=dict(s)
for line in lesi[1:]:
  wrd=line.split()
  k=wrd[1] + "."+wrd[3]
  scoreDict[k]=(scoreDict[k][0],scoreDict[k][1],wrd[15],wrd[13])
scores=scoreDict.values()
scores.sort(reverse=True)
#labels = [sc[1] for sc in scores]
targets,decoys=0.0,0.0
for sc,l,pep,prot in scores:
  if l:
    targets+=1.0
    print "%s\t%f\t%f\t%s\t%s" % ("t_"+str(targets),sc,0.9*decoys/targets,pep,prot)
  else:
    decoys+=1
