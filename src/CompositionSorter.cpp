#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <algorithm>
#include "Globals.h"
#include "Scores.h"
#include "CompositionSorter.h"



// helper for sorting ScoreHolder references by scan and score, to help TDC
bool compareByScanThenScore(const ScoreHolder* lhs, const ScoreHolder* rhs) {
    // Primary sort by scan in ascending order
    if (lhs->pPSM->scan != rhs->pPSM->scan) {
        return lhs->pPSM->scan < rhs->pPSM->scan;
    }
    // Secondary sort by charge in ascending order
    //if (lhs->pPSM->charge != rhs->pPSM->charge) {
    //    return lhs->pPSM->charge < rhs->pPSM->charge;
    //}
    // Thirdly sort by score in descending order
    return lhs->score > rhs->score;
}

void targetDecoyCompetition(std::vector<ScoreHolder*>& scoreHolders) {
    std::sort(scoreHolders.begin(), scoreHolders.end(), compareByScanThenScore);
    cerr << "Before TDC there are " << scoreHolders.size() << " PSMs." << endl;
    // First, sort the vector using the custom comparator

    // Use a map to track scan numbers and keep only the best score for each scan
    std::map<int, ScoreHolder*> bestScoreByScan;
    for (ScoreHolder* holder : scoreHolders) {
        int scan = holder->pPSM->scan;
        if (bestScoreByScan.find(scan) == bestScoreByScan.end() || bestScoreByScan[scan]->score < holder->score) {
            bestScoreByScan[scan] = holder;
        }
    }

    // Clear the original scoreHolders and repopulate it with the best scores from each scan
    scoreHolders.clear();
    for (const auto& pair : bestScoreByScan) {
        scoreHolders.push_back(pair.second);
    }
    cerr << "After TDC there are " << scoreHolders.size() << " PSMs." << endl;
}


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

int CompositionSorter::addPSMs(Scores& scores, bool useTDC) {
    if (useTDC) {
        std::vector<ScoreHolder*> scoreHolders;
        for (auto& scr : scores) {
            scoreHolders.push_back(&scr);
        }
        targetDecoyCompetition(scoreHolders);
        for (auto& pScr : scoreHolders) {
            std::string peptide = pScr->getPSM()->getPeptideSequence();
            std::string signature = generateCompositionSignature(peptide);
            // cerr << "Will push back " << peptide << " with signature " << signature << " there are now " << compositionToPeptidesToScore_[signature][peptide].size() << " such peptides" << endl;
            compositionToPeptidesToScore_[signature][peptide].push_back(pScr);
            // cerr << "Pushed back " << peptide << " with signature " << signature << " there are now " << compositionToPeptidesToScore_[signature][peptide].size() << " such peptides" << endl;
        }
    } else {
        for (auto& scr : scores) {
            std::string peptide = scr.getPSM()->getPeptideSequence();
            std::string signature = generateCompositionSignature(peptide);
            compositionToPeptidesToScore_[signature][peptide].push_back(&scr);
        }
    }
    return 0;
}

int CompositionSorter::sortScorePerPeptide() {
    // Comparator for sorting ScoreHolder references
    auto compareScoreHolder = [](const ScoreHolder* lhs, const ScoreHolder* rhs) {
        return lhs->score > rhs->score; // For descending order
    };

    // Iterating and sorting
    for (auto& [composition, peptideMap] : compositionToPeptidesToScore_) {
        for (auto& [peptide, scores] : peptideMap) {
            std::sort(scores.begin(), scores.end(), compareScoreHolder);
        }
    }
    return 0;
}



int CompositionSorter::inCompositionCompetition(std::vector<ScoreHolder*>& bestScoreHolders, unsigned int decoysPerTarget) {
    // Sort the PSMs for each peptide in decending order of score
    CompositionSorter::sortScorePerPeptide();

    cerr << "inCompositionCompetition, we start with "  << compositionToPeptidesToScore_.size() << " compositions." << endl;

    for (auto& [composition, peptideMap] : compositionToPeptidesToScore_) {
        // cerr << "Composition " << composition << " contains " << peptideMap.size() << " peptides." << endl;
        std::vector<std::vector<ScoreHolder*>> compositionGroups;
        std::vector<ScoreHolder*> targets;
        std::vector<ScoreHolder*> decoys;
        // Add the target peptides
        for (auto& [peptide, scoreHolders] : peptideMap) {
            // cerr << "Peptide, " << peptide << ", with label=" << scoreHolders.front()->label << " has " <<  scoreHolders.size() << " instances." << endl;
            if (scoreHolders.empty()) {
                continue;
            }
            ScoreHolder* firstScoreHolder = scoreHolders.front();
            if (firstScoreHolder->label > 0) {
                targets.push_back(firstScoreHolder);
            } else {
                decoys.push_back(firstScoreHolder);
            }
        }
        // Format tuples of targets and decoys. Now we are not checking if there are enough decoys for each target.
        // TODO: Check if there are enough decoys for each target
        std::reverse(decoys.begin(), decoys.end());
        for (ScoreHolder* target : targets) {
            std::vector<ScoreHolder*> group;
            group.push_back(target); // Add the target ScoreHolder to the group

            // Add decoysPerTarget number of decoys to the group
            for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty(); ++i) {
                group.push_back(decoys.back()); // Add the last decoy (the one first added to the vector)
                decoys.pop_back(); // Remove the added decoy from the decoys list
            }
            compositionGroups.push_back(group); // Add the newly formed group to compositionGroups
        }
        // Take care of the remaining decoys
        for (ScoreHolder* decoy : decoys) {
            std::vector<ScoreHolder*> group;
            group.push_back(decoy); // Add the target ScoreHolder to the group

            // Add decoysPerTarget number of decoys to the group
            for (unsigned int i = 1; i < decoysPerTarget && !decoys.empty(); ++i) {
                group.push_back(decoys.back()); // Add the last decoy (the one first added to the vector)
                decoys.pop_back(); // Remove the added decoy from the decoys list
            }
            compositionGroups.push_back(group); // Add the newly formed group to compositionGroups
        }       
        // Sort each group in compositionGroups based on the ScoreHolder's score
        // and add it to bestScoreHolders
        for (auto& group : compositionGroups) {
            if (group.empty())
                continue; 
            std::sort(group.begin(), group.end(), [](const ScoreHolder* a, const ScoreHolder* b) -> bool {
                return a->score > b->score; // Sort in descending order of score
            });
            bestScoreHolders.push_back(group.front()); // Add the highest-scoring ScoreHolder
            // cerr << "Added pepide, " << group.front()->pPSM->getPeptideSequence() << ", with label=" << group.front()->label << " from a composition group with " << group.size() << " members." << endl; 
        }
    }
    int numTargets = 0;
    for (const auto pSH : bestScoreHolders) {
        if (pSH->label>0) numTargets++;
    }
    if (VERB>1)
        cerr << "inCompositionCompetition, ends with " << bestScoreHolders.size() << " peptides, whereof " << numTargets << " is target peptides." << endl;
    return 0;
}


int CompositionSorter::psmAndPeptide(Scores& scores, std::vector<ScoreHolder *>& winnerPeptides, unsigned int decoysPerTarget) {
    //    psmLevelCompetition(); // Should be handled by the reader?

    // Populate the structure with the winning PSMs
    CompositionSorter::addPSMs(scores);

    // Split out tuples of peptides of identical composition, and select the most high scoring peptide in each tuple
    inCompositionCompetition(winnerPeptides, decoysPerTarget);

    return 0;
}
