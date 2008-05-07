/*******************************************************************************
 Copyright (c) 2008 Oliver Serang

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software. 
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
 
 $Id: ArrayTypes.cpp,v 1.2 2008/05/07 21:25:07 lukall Exp $
 
 *******************************************************************************/

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
