/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

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
#ifndef UNINORMALIZER_H_
#define UNINORMALIZER_H_

class UniNormalizer : public Normalizer { // virtual Normalizer
 public:
  UniNormalizer();
  virtual ~UniNormalizer();
  virtual void setSet(vector<double*> & featuresV,
                      vector<double*> & rtFeaturesV, size_t numFeatures,
                      size_t numRetentionFeatures);
  virtual void updateSet(vector<double*> & featuresV, size_t offset,
                         size_t numFeatures);
  void unnormalizeweight(const vector<double>& in, vector<double>& out);
  void normalizeweight(const vector<double>& in, vector<double>& out);
};

#endif /*UNINORMALIZER_H_*/
