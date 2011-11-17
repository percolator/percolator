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
    cout.precision(10);

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
      //	  cout << "\tinsert " << falsePositiveNames[k];
    }

    //      cout << "Found TP, FP names: " << endl;
    //            cout << truePositiveNames << endl;
    //      cout << falsePositiveNames << endl << endl;

    Array<string> protsAtThreshold;
    string line;

    double prob, lastProb=-1;

    Array<double> empericalList, estimatedList;
    int fpCount = 0, tpCount = 0;

    int numScored = 0;
    Array<string> observedProteins;

    double estFDR = 0.0;
    double empericalFDR = 0.0;
    double totalFDR = 0.0;

    bool scheduledUpdate = false;

    ifstream fin2(argv[2]);

    //while ( cin >> prob && getline(cin, line) )
    while ( fin2 >> prob && getline(fin2, line) )
    {

      //	  cerr << "prob = " << prob << ", \tline = " << line << endl;

      istringstream ist(line);
      ist >> protsAtThreshold;

      //	  cerr << "\tRead " << protsAtThreshold << endl;

      numScored += protsAtThreshold.size();
      observedProteins.append( protsAtThreshold );

      int fpChange = matchCount(falsePosSet, protsAtThreshold);
      int tpChange = matchCount(truePosSet, protsAtThreshold);

      //	  cout << "This iter has fpChange, tpChange = " << fpChange << ", " << tpChange << endl;
      //	  cout << "\tNumber fps, tps = " << fpChange << ", " << tpChange << endl;

      if ( prob != lastProb && lastProb != -1 )
      {
        scheduledUpdate = true;
        //	      cout << "prob has changed from " << lastProb << " to " << prob << endl;
      }

      if ( scheduledUpdate )
      {
        if ( fpChange > 0 || tpChange > 0)
        {
          estimatedList.add(estFDR);
          empericalList.add(empericalFDR);


          //		  cout << "totalFDR = " << totalFDR << endl;
          //		  cout << estFDR << ", " << empericalFDR << "  " << fpCount << " " << tpCount << "   " << prob << endl;

          //		  if ( fpChange > 0 )
          //		    cout << fpCount << " " << tpCount << " " << endl;

          //		  cout << endl;
          //		  cout << "Adding cumulative totals " << fpCount << ", " << tpCount << endl;
          scheduledUpdate = false;

          //		  if ( fpCount > 100 )
          //		    break;

          if ( tpChange > 0 )
          {
            //		      		      cerr << "Adding new TPs at threshold " << prob << endl;
            //		      		      cerr << "\tto make data point with FP, TP = " << fpCount+fpChange << ", " << tpCount+tpChange << endl;
            //		      		      cerr << "\tTPs added were " << matches(truePosSet, protsAtThreshold ) << endl;
          }
        }
      }

      fpCount += fpChange;
      tpCount += tpChange;

      // add 1-prob for all of the non-ignored proteins
      totalFDR += (1-prob) * (fpChange + tpChange);

      estFDR = totalFDR / (fpCount + tpCount);

      empericalFDR = double(fpCount) / (fpCount + tpCount);


      /***
	  if ( fpCount == 10 )
	    {
	      Array<string> tpsAt10 = matches(truePosSet, observedProteins);
	      cout << "AT10: " << tpsAt10 << endl;
	      return 0;
	    }
       ***/


      lastProb = prob;
    }

    //      cerr << "Fin loop" << endl;

    lastProb = prob;

    //      if ( scheduledUpdate )
    {
      estimatedList.add(estFDR);
      empericalList.add(empericalFDR);
    }

    //		cerr << "Final counts are " << fpCount << ", " << tpCount << endl;
    //		cerr << "Overall: " << matchCount(falsePosSet, observedProteins) << ", " << matchCount(truePosSet, observedProteins) << endl;
    //		cerr << "Total number of proteins scored is " << numScored << endl;

    cout << estimatedList << " " << empericalList << endl;
    ofstream fout("/tmp/rlistFDROut.txt");
    fout << estimatedList << endl << empericalList << endl;

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

