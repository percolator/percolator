#!/usr/bin/env python
import sys
ffeatures = open(sys.argv[1],"r")
lf=ffeatures.readlines()
ffeatures.close()
scores=[(float(line.split()[1]),line.split()[0].find("shuffled")<0) for line in lf]
scores.sort(reverse=True)
#labels = [sc[1] for sc in scores]
f,n=0.0,0.0
for sc,l in scores:
  if l:
#    f+=1.0
#    print "1\t%f\t%f" % (sc,n/f)
    print "1"
  else:
#    n+=1.0
#    print "-1\t%f\t%f" % (sc,n/f)
    print "-1"
