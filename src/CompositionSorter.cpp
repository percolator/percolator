#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>
#include "Scores.h"
#include "CompositionSorter.h"


std::string CompositionSorter::generateCompositionSignature(const std::string& peptide) {
    std::map<std::string, int> counts;

    for (size_t i = 0; i < peptide.size(); ++i) {
        if (peptide[i] == '[') {
            size_t j = i + 1;
            while (j < peptide.size() && peptide[j] != ']') {
                ++j;
            }
            // Extract the modification
            counts[peptide.substr(i, j - i + 1)]++;
            // Skip the modification
            i = j;
        } else {
            counts[std::string(1, peptide[i])]++;
        }
    }

    std::string signature;
    for (const auto& [aminoAcid, count] : counts) {
        signature += aminoAcid;
        signature += std::to_string(count);
    }

    return signature;
}

int CompositionSorter::addPSMs(std::vector<ScoreBase>& scores) {
    for (const auto& scr : scores) {
        std::string peptide = scr.getPSM()->getPeptide();
        std::string signature = generateCompositionSignature(peptide);
        if (scr.getPSM()->isTarget()) {
            targetCompositionToPeptides_[signature][peptide].push_back(scr);
        } else {
            decoyCompositionToPeptides_[signature][peptide].push_back(scr);
        }
    }
}


int CompositionSorter::fillPeptidepairs(Scores& scores) {
    // Iterate through the first map
    for (const auto& [composition, targetMap] : targetCompositionToPeptides_) {
        
        // Check if the key exists in the second map
        if (decoyCompositionToPeptides_.count(composition) > 0) {
            const auto& decoyMap = decoyCompositionToPeptides_[composition];

            auto targetIt = targetMap.begin();
            auto decoyIt = decoyMap.begin();

            // Synchronized iteration over inner maps
            while (targetIt != targetInnerMap.end() && decoyIt != decoyInnerMap.end()) {
                auto PeptidepairScore& ppscore(targetIt->second, decoyIt->second, decoyIt->second);
                scores.addScoreHolder(ppscore);
                ++targetIt; ++decoyIt;
            }
            

            // Iterate through the inner map of the first map
            for (const auto& [targetPeptide, targetPSMs] : targetMap) {

                // Check if the inner key exists in the inner map of the second map
                if (decoyInnerMap.count(innerKey) > 0) {
                    const auto& targetVector = targetInnerEntry.second;
                    const auto& decoyVector = decoyInnerMap.at(innerKey);

                    // Now compare vectors
                    for (size_t i = 0; i < targetVector.size() && i < decoyVector.size(); ++i) {
                        // Compare your PSMDescription pointers/objects here
                        if (*targetVector[i] == *decoyVector[i]) {
                            // This indicates that the PSMDescription pointers/objects are equal, you can process them as needed
                            // ... processing code here ...
                            std::cout << "Matching PSMDescription found!" << std::endl;
                        }
                    }
                }
            }
        }
    }

    return 0;
}





    for (const auto& [comp, mp_targ] : targetCompositionToPeptides_) {
        if decoyCompositionToPeptides_.contains(comp) {
            const auto& mp_dec = decoyCompositionToPeptides_[comp]
            // entry.first composition
            // entry.second vector<PSMDescription *>
            for (const auto& [peptide, psmVec] : mp_targ) {
                std::cout << pep << "\n";
            }
            std::cout << "-----\n";
        }
    }

    return 0;
}
