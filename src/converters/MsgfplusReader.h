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

#ifndef MSGFPLUSREADER_H
#define MSGFPLUSREADER_H

#include "MzidentmlReader.h"




class MsgfplusReader : public MzidentmlReader
{
  public:
    MsgfplusReader(ParseOptions *po);

    virtual ~MsgfplusReader();
    bool checkValidity(const std::string &file);
    virtual void searchEngineSpecificParsing(const ::mzIdentML_ns::SpectrumIdentificationItemType & item, int itemCount);
    void addFeatureDescriptions(bool doEnzyme);
    void createPSM(const ::mzIdentML_ns::SpectrumIdentificationItemType & item,
		  ::percolatorInNs::fragSpectrumScan::experimentalMass_type experimentalMass,
		   bool isDecoy, unsigned useScanNumber, boost::shared_ptr<FragSpectrumScanDatabase> database,
		   const std::string & fn);
    double rescaleFragmentFeature(double featureValue, int NumMatchedMainIons);

  protected :
	bool useFragmentSpectrumFeatures;  // Whether to use MS-GF+ high resolution fragmentation spectrum features, default = false
	bool additionalMsgfFeatures;  // Whether the additional features are present in mzid-file, default = false
	int numMatchedIonLimit;  // The number of matched ions required for accurate fragment feature calculation (default 7)
    static const std::map<string,int> msgfplusFeatures; //aux container to map feature name to index
    static const std::map<string,double> msgfplusFeaturesDefaultValue; //aux container to map feature index to feature default value
    const double neutron;  //The difference between C12 and C13
};

#endif // MSGFPLUSREADER_H
