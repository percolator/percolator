// Written by Oliver Serang 2009
// see license for more information

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "GroupPowerBigraph.h"

using namespace std;

int main(int argc, char**argv)
{
  if ( argc == 5 || argc == 6 )
    {
      srand(time(NULL));
      cout.precision(8);

      double gamma = atof(argv[2]);
      double alpha = atof(argv[3]);
      double beta = atof(argv[4]);

      if ( argc == 6 )
	GroupPowerBigraph::LOG_MAX_ALLOWED_CONFIGURATIONS = atof(argv[5]);

      //GroupPowerBigraph gpb( RealRange(alpha, 1, alpha), RealRange(beta, 1, beta), gamma );
      GroupPowerBigraph gpb(argv[1],alpha, beta, gamma );
      //ifstream fin(argv[1]);
      //fin >> gpb;

      gpb.getProteinProbs();
      gpb.printProteinWeights();
    }
  else
    {
      cerr << "usage: Fido graph_file gamma alpha beta [log2_allowed_number_of_connected_states]" << endl;
    }
  return 0;
}

