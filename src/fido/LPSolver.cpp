// Written by Oliver Serang 2009
// see license for more information

#include "LPSolver.h"

void LPSolver::loadArray(istream & fin)
{
  fin >> n;
  fin >> names;

  Array<double> xArray, bArray, fArray;
  Array<Array<double> > AArray;
  fin >> xArray >> fArray >> bArray >> AArray;

  x = xArray;
  f = fArray;
  b = bArray;
  A = AArray;

  iterations = 0;

  addPositivityConstraints();
  normalizeLP();

  //  checkFeasibility();
}

void LPSolver::load(istream & fin)
{
  fin >> n;
  fin >> names;

  fin >> x;
  fin >> f;
  fin >> b;
  fin >> A;

  iterations = 0;

  //  addPositivityConstraints();
  normalizeLP();

  //  checkFeasibility();
}

void LPSolver::normalizeLP()
{
  Matrix newA;
  Vector newB;

  for (int k=0; k<A.numRows(); k++)
    {
      double rowNorm = A[k].norm();

      if ( Vector::sparseChecker.isZero(rowNorm) )
	{
	  cerr << "Error: a row of A has a norm very close to zero" << endl;
	  throw ZeroVectorException();
	}

      double newBValue = b[k] / rowNorm;
      Vector newARow = A[k] / rowNorm;

      newA.add(newARow);
      newB.add(newBValue);
    }

  f = f.normalized();
  A = newA;
  b = newB;
}

void LPSolver::addPositivityConstraints()
{
  Matrix newA;
  Vector newB;

  int k;
  for (k=0; k<n; k++)
    {
      Array<double> slack(n, 0.0);
      slack[k] = 1.0;
      
      newA.add( slack );
      newB.add(0.0);
    }

  for (k=0; k<A.numRows(); k++)
    {
      newA.add( A[k] );
      newB.add( b[k] );
    }

  A = newA;
  b = newB;
}

void LPSolver::printResult() const
{
  cout << "x = " << x << endl;
  cout << "Objective = " << f * x << endl;
  cout << "Iterations = " << iterations << endl;
  cout << endl;
}

void LPSolver::outputMathematica(char * fname) const
{
  ofstream outA(fname);

  outA << "A = " << A.unpack() << ";" << endl;
  outA << "b = " << b.unpack() << ";" << endl;
  outA << "f = " << f.unpack() << ";" << endl;
  outA << "x = " << x.unpack() << ";" << endl;

  outA.close();

  string instr = string("sed s/e/*10^/g ")+string(fname)+string(" > ")+string("/tmp/tmp_t");
  if(system( instr.c_str() ) != 0)
  {
    cerr << "Error: doing system call : " << instr << endl;
  }

  string instr2 = "mv /tmp/tmp_t " + string(fname);
  if(system( instr2.c_str() ) != 0)
  {
    cerr << "Error: doing system call : " << instr2 << endl;
  }


  ofstream outB(fname, ios::app);
  outB << "Timing[LinearProgramming[-f, A, b, Method->\"Simplex\"]]" << endl;
  outB << "f.Last[%]" << endl;
}