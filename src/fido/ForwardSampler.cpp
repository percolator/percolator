// Written by Oliver Serang 2009
// see license for more information

#include "ForwardSampler.h"

ForwardSampler::ForwardSampler():
  hitChecker(1e-10)
{}

ForwardSampler::ForwardSampler(char * fname):
  hitChecker(1e-10)
{
  ifstream fin(fname);
  load(fin);

  addPositivityConstraints();
}

void ForwardSampler::saveArray(ostream & fout)
{
  fout << n << endl;
  fout << names << endl;
  fout << x.unpack() << endl;
  fout << f.unpack() << endl;
  fout << b.unpack() << endl;
  fout << A.unpack() << endl;
}

void ForwardSampler::checkFeasibility() const
{
  Set errors = ((A*x-b) < 0.0);

  if ( ! errors.isEmpty() )
    {
      cerr << "Warning-- didn't start feasible." << endl;
      cerr << "Error constraints are " << errors << endl;
      cerr << "Amount off is " << (A*x-b)[errors] << endl;
      exit(1);
    }  
}

Vector ForwardSampler::solve()
{
  for (;;)
    {
      iterations++;
      descentStep();

      if ( ! randomStep() )
	break;
    }

  return x;
}

Vector ForwardSampler::getProjection() const
{
  // Linear projection based method
  // multiply vector first to ensure O(n^2) for multiplication
  //  return f - S.transpose()*(SSTransInverse*(S*f));

  // Gram-Schmidt based method
  return gs.projectOrthogonal(f);
}

void ForwardSampler::addConstraint(int index)
{
  // Linear projection based method
  //  addConstraintLinearProjection(index);

  // Gram-Schmidt based method
  addConstraintGramSchmidt(index);
}

void ForwardSampler::descentStep()
{
  Vector p = getProjection();

  int count = 0;

  while ( hitChecker.isNonzero( p.norm() / p.size() ) )
    {
      int i = advance(p);

      addConstraint( i );
      p = getProjection();

      count++;
    }
}

Matrix ForwardSampler::getSpanningVectors() const
{
  Matrix sTrans = S.transpose();
  Matrix result =  sTrans.RRE();

  return result;
}

bool ForwardSampler::randomStep()
{
  Set adj = A*x-b == 0;


  Matrix adjA = A[adj];

  Matrix sp = getSpanningVectors();

  // test for optimality--
  // if all spanning vectors have a negative product with f, then
  // you are optimal.
  Set supporting = sp*f > 0;
  Set opposing = sp*f <= 0;

  if ( supporting.isEmpty() )
    {
      // there are only opposing
      // --> optimal; quit
     
      //      cerr << x << endl;
 
      clear();
      return false;
    }

  // now only sample from the ones that can extend

  Set extendable = (adjA * sp.transpose()).transpose().succEq(0.0);

  Set suppExtendable = supporting & extendable;
  Set oppExtendable = opposing & extendable;

  //  Set suppExtendable = supporting;
  //  Set oppExtendable = opposing;

  Array<double> suppCoefs( suppExtendable.size() );
  Array<double> oppCoefs( oppExtendable.size() );

  Random::fillRandomUniform(suppCoefs, 0.0, 1.0);
  Random::fillRandomUniform(oppCoefs, 0.0, 1.0);

  //  cout << "Supporting vector coefficients " << suppCoefs << endl;
  
  Vector direction;
  Vector suppNet, oppNet;

  if ( ! suppExtendable.isEmpty() && ! oppExtendable.isEmpty() )
    {
      suppNet = Vector(suppCoefs)*sp[suppExtendable];
      oppNet = Vector(oppCoefs)*sp[oppExtendable];

      // choose a theta such that ( theta suppNet + (1-theta) oppNet ) * f > 0
      double thetaMin = -oppNet*f/( (suppNet-oppNet)*f );
      double theta = Random::uniform( thetaMin, 1.0 );

      direction = theta * suppNet + (1-theta) * oppNet;

      clear();
    }
  else if ( ! suppExtendable.isEmpty() )
    {
      direction = suppNet;
      clear();
    }
  else
    {
      // performance note: it will be possible to get the appropriate
      // edge in n^2 using matrix inverse updating

      // Bland's rule:
      // I think this should choose the first positive spanning
      // vector, and so should be consistent with Bland's rule
      direction = sp[ supporting[0] ];

      // this should be, by definition, a singleton set
      int dropped = currentGliders[ ( A[currentGliders] * direction > 0 )][0];

      // remove the dropped constraint it from the new glider set
      Set newGliders = currentGliders.without( Set::SingletonSet(dropped) );

      clear();
      currentGliders = newGliders;

      // add everything but the dropped constraint 
      for ( Set::Iterator iter = newGliders.begin(); iter != newGliders.end(); iter++ )
	{
	  addConstraint( *iter );
	}
    }

  int i = advance(direction);

  addConstraint(i);

  return true;
}

void ForwardSampler::clear()
{
  gs = GramSchmidt();
  S = Matrix();
  SSTransInverse = Matrix();
  currentGliders = Set();
}

int ForwardSampler::advance(const Vector & p)
{
  // advances along direction and return the index of the
  // constraint that it hits

  //  cout << "\tgetting opposing..." << endl;
  Set opposing = (A*p) < 0;

  //  cout << "Opposing: " << opposing << endl;
  //  cout << "With A*p: " << A[opposing]*p << endl;

  if ( opposing.isEmpty() )
    {
      // unbounded!
      cerr << "Unbounded LP" << endl;
      exit(1);
    }
  
  // on opposing set:
  // A (x + gamma p) = b
  // find the smallest gamma
  Matrix temp = A[opposing];
  Vector gamma = ( Vector( b[opposing] ) - temp*x ) / ( temp*p );

  double gammaStar = gamma.min();
  Set indices = gamma.argmin();

  x = x + gammaStar * p;

  int firstElement = (opposing[ indices ])[0];

  if ( ! (Set::SingletonSet(firstElement) & currentGliders ).isEmpty() )
    {
      cerr << "\tWarning: adding current glider to set..." << endl;
      cerr << "\t" << firstElement << " added with " << currentGliders.size() << " current members" << endl;
      cerr << "\tIt opposed with " << p*A[firstElement] << " and was reached with gamma = " << gammaStar << endl;
      getchar();
    }

  currentGliders |= Set::SingletonSet( firstElement );

  return firstElement;
}

void ForwardSampler::addConstraintGramSchmidt(int index)
{
  const Vector & row = A[ index ];
  gs.add( row );
  S.add( row );
}

void ForwardSampler::addConstraintLinearProjection(int index)
{
  /***
   // lame method
  S.add(A[index]);
  SSTransInverse = (S * S.transpose() ).inverse();
  return f - S.transpose()*(SSTransInverse*(S*f));
  ***/

  // sherman-morrison method
  const Vector & rowL = A[index];
  S.add( rowL );

  if ( S.numRows() == 1 )
    {
      SSTransInverse = ( S*S.transpose() ).inverse();
    }
  else
    {
      // compute the updated inverse of (S * S^T)
      SSTransInverse.embedInIdentity();

      Vector uA( SSTransInverse.numRows() );
      uA.addElement( SSTransInverse.numRows()-1 , 1.0 );

      Vector uB = S*rowL;
      
      SSTransInverse = ShermanMorrison(SSTransInverse, uB, uA);
      SSTransInverse = ShermanMorrison(SSTransInverse, uA, uB);

      // remove the corner element, which is added twice

      SSTransInverse = ShermanMorrison(SSTransInverse, uA, (-1.0 - 1.0*rowL*rowL)*uA);
    }
}
