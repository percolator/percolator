Array<int> seq(int lowest, int highest)
{
  Array<int> result( highest - lowest + 1 );
  
  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + k;
    }

  return result;
}

Vec seq(double lowest, double highest, double step)
{
  Array<double> result( int( (highest - lowest) / step ) + 1 );

  int k;
  for (k=0; k<result.size(); k++)
    {
      result[k] = lowest + step * k;
    }

  return result;
}
