// Written by Oliver Serang 2009
// see license for more information

#ifndef _PivdoSplitter_H
#define _PivdoSplitter_H

#include "BasicBigraph.h"

class PivdoSplitter : public BasicBigraph
{
 public:
  void outputPivdo(ostream & os) const;
  
  void outputSplitPivdo(int iter, int numberSplits, char * destPath);
  
  PivdoSplitter() {}
 PivdoSplitter(const BasicBigraph & bb): BasicBigraph(bb)
    {
    }
};

#endif
