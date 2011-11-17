// Written by Oliver Serang 2009
// see license for more information

#include "PivdoSplitter.h"

void PivdoSplitter::outputPivdo(ostream & os) const
{
  for (int k=0; k<PSMsToProteins.size(); k++)
    {
      os << "e " << PSMsToProteins.names[k] << endl;

      const Set & s = PSMsToProteins.associations[k];
      for (Set::Iterator iter = s.begin(); iter != s.end(); iter++)
	{
	  os << "r " << proteinsToPSMs.names[ *iter ] << endl;
	}

      os << "p " << PSMsToProteins.weights[k] << endl;
    }
}

void PivdoSplitter::outputSplitPivdo(int fold, int numFolds, char* destPath) 
{
  // note: shuffling these later may be a good idea


  // this is important-- it guarantees shuffling will not make the validation sections overlap
  srand(0);

  PeptideThreshold = -1;

  prune();

  Array<BasicBigraph> parts = partitionSections();
  Array<PivdoSplitter> psParts = Array<PivdoSplitter>(parts.size());

  // shuffle these subgraphs
  Array<int> shuffledIndices(parts.size());
  
  int i;
  for (i=0; i<psParts.size(); i++)
    {
      shuffledIndices[i] = i;
    }
  for (i=0; i<psParts.size(); i++)
    {
      int otherInd = rand() % psParts.size();

      swap( shuffledIndices[i], shuffledIndices[otherInd] );
    }

  for (i=0; i<psParts.size(); i++)
    psParts[i] = PivdoSplitter(parts[ shuffledIndices[i] ]);

  // split for the fold x-val

  ostringstream trainName;
  trainName << destPath << "train_" << fold << ".pivdo";
  ofstream trainFout(trainName.str().c_str());

  ostringstream testName;
  testName << destPath << "test_" << fold << ".pivdo";
  ofstream testFout(testName.str().c_str());

  int numProtsSeen = 0;
  int currentFold = 0;

  for ( int k=0; k<parts.size(); k++)
    {
      //      cout << currentFold << " " << fold << endl;

      //      cout << "numProtsSeen: " << numProtsSeen << "  -->  " << numProtsSeen+parts[k].proteinsToPSMs.size() << endl;

      if ( currentFold == fold )
	psParts[k].outputPivdo(testFout);
      else
	psParts[k].outputPivdo(trainFout);

      numProtsSeen += parts[k].proteinsToPSMs.size();
      if ( numProtsSeen >= ( ((currentFold+1) * proteinsToPSMs.size()) / numFolds ) )
	{
	  //	  cout << "Changing fold on this one.." << endl << endl;
	  currentFold++;
	}
    }
}


