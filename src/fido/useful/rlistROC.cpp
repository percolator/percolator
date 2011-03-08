#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"
#include <set>
//#include "Hash_Set.h"



using namespace std;

#include <ext/hash_set>

using namespace std;

namespace __gnu_cxx
{
template<> struct hash< std::string >
{
    size_t operator()( const std::string & x ) const
    {
      return hash< const char* >()( x.c_str() );
    }
};
}


int matchCount( const __gnu_cxx::hash_set<string> & positiveNames, const Array<string> & atThreshold )
{
  int count = 0;

  for (int k=0; k<atThreshold.size(); k++)
  {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
      count++;
  }

  return count;
}

Array<string> matches( const __gnu_cxx::hash_set<string> & positiveNames, const Array<string> & atThreshold )
    {
  Array<string> result;
  for (int k=0; k<atThreshold.size(); k++)
  {
    if ( positiveNames.count( atThreshold[k] ) > 0 )
      result.add( atThreshold[k] );
  }
  return result;
    }

int main(int argc, char**argv)
{
  if ( argc == 3 )
  {
    ifstream fin(argv[1]);

    Array<string> truePositiveNames, falsePositiveNames;
    fin >> truePositiveNames;
    fin >> falsePositiveNames;

    __gnu_cxx::hash_set<string> truePosSet(truePositiveNames.size()), falsePosSet(falsePositiveNames.size());

    int k;
    for (k=0; k<truePositiveNames.size(); k++)
    {
      truePosSet.insert( truePositiveNames[k] );
      //	  cout << "\tinsert " << truePositiveNames[k];
    }
    for (k=0; k<falsePositiveNames.size(); k++)
    {
      falsePosSet.insert( falsePositiveNames[k] );
    }

    Array<string> protsAtThreshold;
    string line;

    double prob, lastProb=-1;

    Array<int> fps, tps;
    int fpCount = 0, tpCount = 0;

    int numScored = 0;
    Array<string> observedProteins;

    fps.add(0);
    tps.add(0);
    bool scheduledUpdate = false;

    double totalFDR = 0.0, estFDR = 0.0;

    ifstream fin2(argv[2]);

    //while ( cin >> prob && getline(cin, line) )
    while ( fin2 >> prob && getline(fin2, line) )
    {
      istringstream ist(line);
      ist >> protsAtThreshold;

      numScored += protsAtThreshold.size();
      observedProteins.append( protsAtThreshold );

      int fpChange = matchCount(falsePosSet, protsAtThreshold);
      int tpChange = matchCount(truePosSet, protsAtThreshold);

      // for different style of grading, recognizing only target
      // proteins in merged protein groups
      /**
	  if ( tpChange > 0 && fpChange > 0 )
	    {
	      cerr << "Note: When FPs has already been at " << fps.back() << "(est FDR = " << estFDR << "), losing " << fpChange << " false positives" << endl;
	      fpChange = 0;
	    }
       **/
      // end

      if ( prob != lastProb && lastProb != -1 )
      {
        scheduledUpdate = true;
      }

      if ( scheduledUpdate )
      {
        fps.add( fpCount );
        tps.add( tpCount );
        scheduledUpdate = false;

        //	      if ( fpCount > 100 )
        //		break;

        totalFDR += (1-prob) * (fpChange + tpChange);
        estFDR = totalFDR / (fpCount + tpCount);

      }

      fpCount += fpChange;
      tpCount += tpChange;


      lastProb = prob;
    }

    lastProb = prob;

    fps.add( fpCount );
    tps.add( tpCount );

    fps.add( falsePosSet.size() );
    tps.add( truePosSet.size() );

    cout << fps << " " << tps << endl;
    ofstream fout("/tmp/fido/rlistROCOut.txt");
    fout << fps << endl << tps << endl;

    // closing file streams
    fout.close();
    fin2.close();
    fin.close();
  }
  else
  {
    cerr << "usage: <True Positives and Overall Array File>" << endl << endl;
  }
  return 0;
}

