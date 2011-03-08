#include "Matrix.h"
#include "Vector.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

double THRESH;

double squareAntiderivativeAt(double m, double b, double xVal)
{
  // turn into ux^2+vx+t
  double u = m*m;
  double v = 2*m*b;
  double t = b*b;

  //  cout << "\t\tsquareAntiderivativeAt " << xVal << " is " << u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal << endl;
  //  cout << "u, v, t = " << u << " " << v << " " << t << endl;

  return u*xVal*xVal*xVal/3.0 + v*xVal*xVal/2.0 + t*xVal;
}

double antiderivativeAt(double m, double b, double xVal)
{
  return m*xVal*xVal/2.0 + b*xVal;
}

double area(double x1, double y1, double x2, double y2)
{
  double m = (y2-y1)/(x2-x1);
  double b = y1-m*x1;

  //  cout << "\tarea from region " << x1 << " " << x2 << "  " << squareAntiderivativeAt(m, b, min(THRESH, x2) ) - squareAntiderivativeAt(m, b, x1) << endl;
  //  cout << "\tusing m = " << y2-y1 << "/" << x2-x1 << " = "  << m << endl;

  return squareAntiderivativeAt(m, b, min(THRESH, x2) ) - squareAntiderivativeAt(m, b, x1);
}


int main(int argc, char**argv)
{
  if ( argc == 3 )
  {
    THRESH = atof(argv[1]);

    string both = argv[2];
    //istringstream istBoth(both);
    ifstream istBoth(argv[2]);

    Array<double> estFDR, empFDR;

    istBoth >> estFDR;
    istBoth >> empFDR;

    if ( estFDR.size() != empFDR.size() )
    {
      cerr << "Error: Arrays do not have the same size" << endl;
      exit(1);
    }

    Vector diff = Vector(estFDR) - Vector(empFDR);

    double tot = 0.0;

    int k;
    for (k=0; k<diff.size()-1; k++)
    {
      // stop if no part of the estFDR is < threshold
      if ( estFDR[k] >= THRESH )
      {
        if ( k == 0 )
          tot = 1.0 / 0.0;

        break;
      }

      tot += area(estFDR[k], diff[k], estFDR[k+1], diff[k+1]);
    }

    double xRange = min(THRESH, estFDR[k]) - estFDR[0];

    if ( isinf(tot) )
      cout << tot << endl;
    else
      cout << tot / xRange;
  }
  else
  {
    cerr << "usage: MSE FDR_threshold \"<array double estimated FDR> <array double emperical FDR>\"" << endl << endl;
  }
  return 0;
}

