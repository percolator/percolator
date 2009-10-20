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
using namespace std;

class Normalizer
{
public:
	virtual ~Normalizer();
    virtual void setSet(vector<double *> & featuresV,vector<double *> & rtFeaturesV, size_t numFeatures, size_t numRetentionFeatures){;}
//    virtual void setPsmSet(vector<PSMDescription> & psms, size_t noFeat){;}
//    void normalizeSet(vector<PSMDescription> & psms);
    void normalizeSet(vector<double *> & featuresV,vector<double *> & rtFeaturesV);
    void normalize(const double * in, double * out, size_t offset, size_t numFeatures);
    double normalize(const double in, size_t index) { return (in-sub[index])/div[index]; }
    virtual void unnormalizeweight(const vector<double>& in,vector<double>& out){;}
    virtual void normalizeweight(const vector<double>& in,vector<double>& out){;}
    static Normalizer * getNormalizer();
    static void setType(int type);
	const static int UNI = 0;
	const static int STDV = 1;
	void resizeVecs(size_t size) {sub.resize(size); div.resize(size); }
	void setNumberRetentionFeatures(size_t numRF) {numRetentionFeatures = numRF; }
	double* getSub() {return &sub[0];}
	double* getDiv() {return &div[0];}
	size_t*  getNumRetFeatures() {return & numRetentionFeatures;}
	void printNumRetFeatures() {cout << "There are " << numRetentionFeatures << endl; }
	void printSub()  {for (int i = 0; i < numRetentionFeatures; i++) cout << sub[i] << " "; cout << endl;}
	void printDiv()  {for (int i = 0; i < numRetentionFeatures; i++) cout << div[i] << " "; cout << endl;}
protected:
    Normalizer();
    static Normalizer * theNormalizer;
	static int subclass_type;
    size_t numFeatures, numRetentionFeatures;
    vector<double> sub;
    vector<double> div;
};

#endif /*NORMALIZER_H_*/
