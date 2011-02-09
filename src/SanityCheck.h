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
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_

class Scores;
class Normalizer;

class SanityCheck {
  public:
    SanityCheck();
    virtual ~SanityCheck();
    static SanityCheck* initialize(string otherCall);
    void readWeights(istream& weightStream, vector<double>& w);
    int getInitDirection(vector<Scores>& testset,
                         vector<Scores> &trainset, Normalizer* pNorm,
                         vector<vector<double> >& w, double test_fdr);
    virtual bool validateDirection(vector<vector<double> >& w);
    void resetDirection(vector<vector<double> >& w);

    static void setInitWeightFN(string fn) {
      initWeightFN = fn;
    }
    static void setInitDefaultDir(int dir) {
      initDefaultDir = dir;
    }
    static void setOverrule(bool orl) {
      overRule = orl;
    }
  protected:
    virtual void getDefaultDirection(vector<vector<double> >& w);
    int initPositives;
    double fdr;
    static bool overRule;
    static string initWeightFN;
    static int initDefaultDir; // Default Direction, 0=do not use,
    // positive integer = feature number,
    // negative integer = lower score better
    vector<Scores> *pTestset, *pTrainset;
};

#endif /*SANITYCHECK_H_*/
