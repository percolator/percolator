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

/*
 * Unit tests for the peptide selection of the RESET algorithm
 * (CompositionSorter): the winners, and their order, are compared with the
 * earlier implementations, which are kept here as references.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <deque>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "CompositionSorter.h"
#include "Globals.h"
#include "PsmStore.h"
#include "Scores.h"

namespace {

// Random PSMs whose peptides repeat, with coarse scores so that some tie.
// The peptides (with flanks, "K.PEPTIDE.R") are kept in a deque, whose
// elements do not move, since the PSMs only view them.
class RandomPsms {
  public:
    RandomPsms(unsigned int seed, int numPsms, int numSequences,
               const std::string& alphabet, bool withMods)
        : seed_(seed) {
        std::vector<std::string> sequences;
        for (int k = 0; k < numSequences; ++k) {
            std::string sequence;
            const int length = 2 + static_cast<int>(next() % 12);
            for (int i = 0; i < length; ++i) {
                sequence += alphabet[next() % alphabet.size()];
                if (withMods && next() % 8 == 0) {
                    static const char* const mods[] = {
                        "[+16]", "[+57.02]", "[-17]", "[]", "[[x]", "[+16]"};
                    sequence += mods[next() % 6];
                }
            }
            if (withMods && next() % 20 == 0) sequence += "[unclosed";
            sequences.push_back(sequence);
        }
        psms_.resize(numPsms);
        for (int i = 0; i < numPsms; ++i) {
            peptides_.push_back("K." + sequences[next() % sequences.size()] +
                                ".R");
            store_.storeStrings(psms_[i], "", peptides_.back(), {});
            psms_[i].scan = static_cast<unsigned int>(next() % (numPsms + 1));
            const bool isTarget = next() % 2 == 0;
            scores_.addScoreHolder(ScoreHolder(
                static_cast<double>(next() % 50) / 4.0,
                isTarget ? LabelType::TARGET : LabelType::DECOY, &psms_[i]));
        }
        scores_.recalculateSizes();
    }
    const Scores& scores() const { return scores_; }

  private:
    unsigned int next() {
        seed_ = seed_ * 1103515245u + 12345u;
        return seed_ >> 8;
    }
    unsigned int seed_;
    std::deque<std::string> peptides_;
    PsmStore store_;  // holds the peptides of psms_
    std::vector<PSMDescription> psms_;
    Scores scores_{false};
};

// CompositionSorter::psmsOnly before it used std::string_view keys.
std::vector<ScoreHolder> referencePsmsOnly(const Scores& scores) {
    std::unordered_map<std::string, ScoreHolder*> bestScoreHolders;
    for (const ScoreHolder& sh : scores) {
        const std::string& peptide = sh.getPSM()->getPeptideSequence();
        if (bestScoreHolders.find(peptide) == bestScoreHolders.end() ||
            sh.score > bestScoreHolders[peptide]->score) {
            bestScoreHolders[peptide] = const_cast<ScoreHolder*>(&sh);
        }
    }
    std::vector<ScoreHolder> winners;
    for (const auto& entry : bestScoreHolders) {
        winners.push_back(*(entry.second));
    }
    return winners;
}

void expectSameWinners(const std::vector<ScoreHolder>& expected,
                       const Scores& winners, const std::string& what) {
    ASSERT_EQ(expected.size(), static_cast<size_t>(winners.end() -
                                                   winners.begin()))
        << what;
    size_t k = 0;
    for (const ScoreHolder& sh : winners) {
        EXPECT_EQ(expected[k].pPSM, sh.pPSM) << what << " winner " << k;
        EXPECT_EQ(expected[k].score, sh.score) << what << " winner " << k;
        EXPECT_EQ(expected[k].label, sh.label) << what << " winner " << k;
        ++k;
    }
}

class CompositionSorterTest : public ::testing::Test {
  protected:
    void SetUp() override {
        origVerbose_ = Globals::getInstance()->getVerbose();
        Globals::getInstance()->setVerbose(0);
    }
    void TearDown() override {
        Globals::getInstance()->setVerbose(origVerbose_);
    }
  private:
    int origVerbose_ = 0;
};

}  // namespace

// psmsOnly selects the same best PSM per peptide, in the same order, as
// the implementation with std::string keys.
TEST_F(CompositionSorterTest, PsmsOnlyKeepsTheWinnersAndTheirOrder)
{
    for (int round = 0; round < 8; ++round) {
        RandomPsms psms(17u + round, 200 + 900 * round, 50 + 300 * round,
                        "ACDEFGHIKLMNPQRSTVWY", round % 2 == 1);
        Scores winners(false);
        CompositionSorter::psmsOnly(psms.scores(), winners);
        expectSameWinners(referencePsmsOnly(psms.scores()), winners,
                          "round " + std::to_string(round));
        EXPECT_EQ(winners.size(),
                  static_cast<unsigned int>(winners.end() - winners.begin()));
    }
}

// defined in CompositionSorter.cpp
void targetDecoyCompetition(std::vector<ScoreHolder*>& scoreHolders);

namespace {

// CompositionSorter::generateCompositionSignature with a std::map per
// peptide, before the counts were kept in a table.
std::string referenceSignature(const std::string& peptide) {
    std::map<std::string, int> counts;
    for (size_t i = 0; i < peptide.size(); ++i) {
        if (peptide[i] == '[') {
            size_t j = i + 1;
            while (j < peptide.size() && peptide[j] != ']') {
                ++j;
            }
            counts[peptide.substr(i, j - i + 1)]++;
            i = j;
        } else {
            counts[std::string(1, peptide[i])]++;
        }
    }
    std::string signature;
    for (const auto& aminoAcidCountPair : counts) {
        signature += aminoAcidCountPair.first;
        signature += std::to_string(aminoAcidCountPair.second);
    }
    return signature;
}

// The composition competition of CompositionSorter before it kept the PSMs
// in one sorted vector: a hash map of compositions, each with a std::map of
// its peptides and their PSMs.
std::vector<ScoreHolder> referenceCompositionWinners(
        const Scores& scores, unsigned int decoysPerTarget, bool useTDC) {
    std::unordered_map<std::string,
                       std::map<std::string, std::vector<ScoreHolder*>>,
                       DJB2Hash> compositionToPeptidesToScore;
    std::vector<ScoreHolder*> scoreHolders;
    for (const auto& scr : scores) {
        scoreHolders.push_back(const_cast<ScoreHolder*>(&scr));
    }
    if (useTDC) targetDecoyCompetition(scoreHolders);
    for (ScoreHolder* pScr : scoreHolders) {
        std::string peptide = pScr->getPSM()->getPeptideSequence();
        std::string signature = referenceSignature(peptide);
        compositionToPeptidesToScore[signature][peptide].push_back(pScr);
    }
    auto compareScoreHolder = [](const ScoreHolder* lhs,
                                 const ScoreHolder* rhs) {
        return *lhs > *rhs;
    };
    for (auto& compositionToPeptides : compositionToPeptidesToScore) {
        for (auto& peptideToScores : compositionToPeptides.second) {
            std::sort(peptideToScores.second.begin(),
                      peptideToScores.second.end(), compareScoreHolder);
        }
    }
    std::vector<ScoreHolder> winners;
    for (auto& compositionToPeptides : compositionToPeptidesToScore) {
        std::vector<std::vector<ScoreHolder*>> compositionGroups;
        std::vector<ScoreHolder*> targets;
        std::vector<ScoreHolder*> decoys;
        for (auto& peptideToScores : compositionToPeptides.second) {
            ScoreHolder* first = peptideToScores.second.front();
            if (first->isTarget()) {
                targets.push_back(first);
            } else {
                decoys.push_back(first);
            }
        }
        std::reverse(decoys.begin(), decoys.end());
        for (ScoreHolder* target : targets) {
            std::vector<ScoreHolder*> group;
            group.push_back(target);
            for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty();
                 ++i) {
                group.push_back(decoys.back());
                decoys.pop_back();
            }
            compositionGroups.push_back(group);
        }
        for (ScoreHolder* decoy : decoys) {
            std::vector<ScoreHolder*> group;
            group.push_back(decoy);
            for (unsigned int i = 1; i < decoysPerTarget && !decoys.empty();
                 ++i) {
                group.push_back(decoys.back());
                decoys.pop_back();
            }
            compositionGroups.push_back(group);
        }
        for (auto& group : compositionGroups) {
            std::sort(group.begin(), group.end(), compareScoreHolder);
            winners.push_back(*(group.front()));
        }
    }
    return winners;
}

}  // namespace

// The signature lists residues and modifications in the order of the
// std::map with std::string keys: also lower-case letters, digits and
// other characters on both sides of '[', bytes above 127, unclosed and
// empty modifications, and '[' inside a modification.
TEST_F(CompositionSorterTest, SignatureMatchesTheMapOfStrings)
{
    CompositionSorter sorter;
    const std::string fixed[] = {
        "", "A", "PEPTIDE", "PEPT[+80]IDE[+80]K", "M[+16]M[+16]M",
        "AC[+57.02]DC[+57.02]", "n[+42]PEPTIDE", "[unclosed", "AB[",
        "[]A[]", "[[x]]", "Z[A]a[z]", "\xc3\xa9t\xc3\xa9", "-.0_9~"};
    for (const std::string& peptide : fixed) {
        EXPECT_EQ(referenceSignature(peptide),
                  sorter.generateCompositionSignature(peptide))
            << "peptide " << peptide;
    }
    unsigned int seed = 31u;
    auto next = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return seed >> 8;
    };
    const std::string alphabet = "ACDEFGHIKLMNPQRSTVWYacm0129-.[]Z^_`{\x7f";
    for (int round = 0; round < 20000; ++round) {
        std::string peptide;
        const int length = static_cast<int>(next() % 25);
        for (int i = 0; i < length; ++i) {
            if (next() % 40 == 0) {
                peptide += static_cast<char>(128 + next() % 128);
            } else {
                peptide += alphabet[next() % alphabet.size()];
            }
        }
        ASSERT_EQ(referenceSignature(peptide),
                  sorter.generateCompositionSignature(peptide))
            << "peptide " << peptide;
    }
}

// psmAndPeptide (--composition-match) selects the same winners, in the
// same order, as the implementation with nested maps, for one and two
// decoys per target, and the same holds with target-decoy competition in
// addPSMs.
TEST_F(CompositionSorterTest, CompositionMatchKeepsTheWinnersAndTheirOrder)
{
    for (int round = 0; round < 8; ++round) {
        // a small alphabet, so that peptides share compositions
        RandomPsms psms(101u + round, 300 + 1000 * round, 80 + 400 * round,
                        round % 2 ? "ACDEK" : "ACDEFGHIKLMNPQRSTVWY",
                        round % 3 != 0);
        for (unsigned int decoysPerTarget : {1u, 2u}) {
            const std::string what = "round " + std::to_string(round) +
                " decoysPerTarget " + std::to_string(decoysPerTarget);
            Scores winners(false);
            CompositionSorter sorter;
            sorter.psmAndPeptide(psms.scores(), winners, decoysPerTarget);
            expectSameWinners(
                referenceCompositionWinners(psms.scores(), decoysPerTarget,
                                            false),
                winners, what);
            EXPECT_EQ(winners.size(), static_cast<unsigned int>(
                                          winners.end() - winners.begin()));
        }
        Scores tdcWinners(false);
        CompositionSorter tdcSorter;
        tdcSorter.addPSMs(psms.scores(), true);
        tdcSorter.inCompositionCompetition(tdcWinners, 1);
        expectSameWinners(referenceCompositionWinners(psms.scores(), 1, true),
                          tdcWinners, "TDC round " + std::to_string(round));
    }
}
