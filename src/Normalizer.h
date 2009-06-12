/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

 *******************************************************************************/
#ifndef NORMALIZER_H_
#define NORMALIZER_H_
class DataSet;
#include <set>
using namespace std;

class Normalizer
{
public:
	virtual ~Normalizer();
    virtual void setSet(set<DataSet *> & setVec, size_t numFeatures, size_t numRetentionFeatures){;}
    void normalizeSet(set<DataSet *> & setVec);
    void normalize(const double * in, double * out, size_t offset, size_t numFeatures);
    double normalize(const double in, size_t index) { return (in-sub[index])/div[index]; }
    virtual void unnormalizeweight(const vector<double>& in,vector<double>& out){;}
    virtual void normalizeweight(const vector<double>& in,vector<double>& out){;}
    static Normalizer * getNormalizer();
    static void setType(int type);
	const static int UNI = 0;
	const static int STDV = 1;
protected:
    Normalizer();
    static Normalizer * theNormalizer;
	static int subclass_type;
    size_t numFeatures, numRetentionFeatures;
    vector<double> sub;
    vector<double> div;
};

#endif /*NORMALIZER_H_*/
