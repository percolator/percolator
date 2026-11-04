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

#include <deque>
#include <string>
#include <unordered_map>
#include <vector>

#include "CompositionSorter.h"
#include "Globals.h"
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
            psms_[i] = PSMDescription(peptides_.back());
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
