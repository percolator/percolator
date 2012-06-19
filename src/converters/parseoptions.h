/*
    Copyright 2012 <copyright holder> <email>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#ifndef PARSEOPTIONS_H
#define PARSEOPTIONS_H

#include <iostream>
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
    reversedFeaturePattern(""){};
    bool calcQuadraticFeatures;
    bool calcAAFrequencies;
    bool calcPTMs;
    bool calcDOC;
    bool isotopeMass;
    int hitsPerSpectrum;
    bool pngasef;
    bool iscombined;
    std::string reversedFeaturePattern;
    std::map<char, int> ptmScheme;
};

#endif // PARSEOPTIONS_H
