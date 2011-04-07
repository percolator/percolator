#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include "Vector.h"

using namespace std;

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2, int N)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  return antiderivativeAt(m, b, min(double(N), x2) ) - antiderivativeAt(m, b, x1);
}

int main(int argc, char**argv)
{
  if ( argc == 3 )
  {
    int N=atoi(argv[1]);

    Array<double> fpArray, tpArray;
    ifstream fin(argv[2]);
    fin >> fpArray >> tpArray;

    double rocN = 0.0;

    if ( fpArray.back() < N )
    {
      cerr << "There are not enough false positives; needed " << N << " and was only given " << fpArray.back() << endl << endl;
      exit(1);
    }

    for (int k=0; k<fpArray.size()-1; k++)
    {
      // find segments where the fp value changes

      if ( fpArray[k] >= N )
        break;

      if ( fpArray[k] != fpArray[k+1] )
      {
        // this line segment is a function

        double currentArea = area(fpArray[k], tpArray[k], fpArray[k+1], tpArray[k+1], N);
        rocN += currentArea;
      }
    }

    double roc50 = rocN / (N * tpArray.back());
    cout << roc50;
    ofstream fout("/tmp/roc50Out.txt");
    fout << roc50;

    // closing file streams
    fout.close();
  }
  else
  {
    cerr << "usage: ROCN n < roc_file" << endl << endl;
  }
  return 0;
}

