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
#ifndef PARSEOPTIONS_H
#define PARSEOPTIONS_H

#include <string>
#include <map>

class ParseOptions
{
  
 public:
    ParseOptions(): calcQuadraticFeatures(false),
    calcAAFrequencies(false),
    calcPTMs(false),
    isotopeMass(false),
    pngasef(false),
    calcDOC(false),
    hitsPerSpectrum(1),
    iscombined(false),
    reversedFeaturePattern("random"),
    peptidelength(6),
    monoisotopic(false),
    boost_serialization(true),
    targetFN(""),
    decoyFN(""),
    call(""),
    spectrumFN(""),
    xmlOutputFN(""),
    minmass(400),
    maxmass(6000),
    maxpeplength(40),
    missed_cleavages(0),
    targetDb(""),
    decoyDb(""),
    readProteins(false)
    {};
    bool calcQuadraticFeatures;
    bool calcAAFrequencies;
    bool calcPTMs;
    bool calcDOC;
    bool isotopeMass;
    int hitsPerSpectrum;
    unsigned peptidelength;
    bool pngasef;
    bool iscombined;
    bool monoisotopic;
    bool boost_serialization;
    std::string reversedFeaturePattern;
    std::string targetFN;
    std::string decoyFN;
    std::string spectrumFN;
    std::string call;
    std::string xmlOutputFN;
    std::map<char, int> ptmScheme;
    double minmass;
    double maxmass;
    unsigned maxpeplength;
    unsigned missed_cleavages;
    std::string targetDb;
    std::string decoyDb;
    bool readProteins;
};

#endif // PARSEOPTIONS_H
