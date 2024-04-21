#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <algorithm>
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

int CompositionSorter::addPSMs(Scores& scores) {
    for (const auto& scr : scores) {
        std::string peptide = scr.getPSM()->getPeptideSequence();
        std::string signature = generateCompositionSignature(peptide);
        compositionToPeptidesToScore_[signature][peptide].push_back(&scr);
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

// Scores
int CompositionSorter::inCompositionCompetition(Scores& bestScoreHolders, unsigned int decoysPerTarget) {

    for (auto& [composition, peptideMap] : compositionToPeptidesToScore_) {
        std::vector<std::vector<const ScoreHolder*>> compositionGroups;
        std::vector<const ScoreHolder*> targets;
        std::vector<const ScoreHolder*> decoys;
        // Add the target peptides
        for (auto& [peptide, scoreHolders] : peptideMap) {
            if (scoreHolders.empty()) {
                continue;
            }
            const ScoreHolder* firstScoreHolder = scoreHolders.front();
            if (firstScoreHolder->label > 0) {
                targets.push_back(firstScoreHolder);
            } else {
                decoys.push_back(firstScoreHolder);
            }
        }
        // Format tuples of targets and decoys. Now we are not checking if there are enough decoys for each target.
        // TODO: Check if there are enough decoys for each target
        for (const ScoreHolder* target : targets) {
            std::vector<const ScoreHolder*> group;
            group.push_back(target); // Add the target ScoreHolder to the group

            // Add decoysPerTarget number of decoys to the group
            for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty(); ++i) {
                group.push_back(decoys.back()); // Add the last decoy to ensure unique selection if decoys are not repeated
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
            bestScoreHolders.addScoreHolder(*(group.front())); // Add the highest-scoring ScoreHolder
        }
     }
    // Select the first ScoreHolder in each group in compositionGroups and feed it to the bestScoreHolders vector
    return 0;
}


int CompositionSorter::inCompositionCompetition(std::vector<ScoreHolder*>& bestScoreHolders, unsigned int decoysPerTarget) {

    for (auto& [composition, peptideMap] : compositionToPeptidesToScore_) {
        std::vector<std::vector<ScoreHolder*>> compositionGroups;
        std::vector<ScoreHolder*> targets;
        std::vector<ScoreHolder*> decoys;
        // Add the target peptides
        for (auto& [peptide, scoreHolders] : peptideMap) {
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
        for (ScoreHolder* target : targets) {
            std::vector<ScoreHolder*> group;
            group.push_back(target); // Add the target ScoreHolder to the group

            // Add decoysPerTarget number of decoys to the group
            for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty(); ++i) {
                group.push_back(decoys.back()); // Add the last decoy to ensure unique selection if decoys are not repeated
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
        }
     }
    // Select the first ScoreHolder in each group in compositionGroups and feed it to the bestScoreHolders vector
    return 0;
}


int CompositionSorter::psmAndPeptide(Scores& scores, Scores& winnerPeptides, unsigned int decoysPerTarget) {
    //    psmLevelCompetition(); // Should be handled by the reader?

    // Populate the structure with the winning PSMs
    CompositionSorter::addPSMs(scores);

    // Sort the PSMs for each peptide in decending order of score
    CompositionSorter::sortScorePerPeptide();

    // Split out tuples of peptides of identical composition, and select the most high scoring peptide in each tuple
    inCompositionCompetition(winnerPeptides, decoysPerTarget);

    return 0;
}

int CompositionSorter::psmAndPeptide(Scores& scores, std::vector<ScoreHolder *>& winnerPeptides, unsigned int decoysPerTarget) {
    //    psmLevelCompetition(); // Should be handled by the reader?

    // Populate the structure with the winning PSMs
    CompositionSorter::addPSMs(scores);

    // Sort the PSMs for each peptide in decending order of score
    CompositionSorter::sortScorePerPeptide();

    // Split out tuples of peptides of identical composition, and select the most high scoring peptide in each tuple
    inCompositionCompetition(winnerPeptides, decoysPerTarget);

    return 0;
}
