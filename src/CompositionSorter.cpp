#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CompositionSorter.h"
#include "Globals.h"
#include "ParallelSort.h"
#include "Scores.h"

// helper for sorting ScoreHolder references by scan and score, to help TDC
bool compareByScanThenScore(const ScoreHolder* lhs, const ScoreHolder* rhs) {
  // sort by file, scan, mass in ascending order
  if (lhs->pPSM->specFileNr != rhs->pPSM->specFileNr) {
    return lhs->pPSM->specFileNr < rhs->pPSM->specFileNr;
  }
  if (lhs->pPSM->scan != rhs->pPSM->scan) {
    return lhs->pPSM->scan < rhs->pPSM->scan;
  }
  if (lhs->pPSM->expMass != rhs->pPSM->expMass) {
    return lhs->pPSM->expMass < rhs->pPSM->expMass;
  }
  // Then sort by score in descending order, ties in the total order of
  // ScoreHolder so that the result does not depend on the sort algorithm
  return *lhs > *rhs;
}

void targetDecoyCompetition(std::vector<ScoreHolder*>& scoreHolders) {
  parallelSort(scoreHolders.begin(), scoreHolders.end(),
               compareByScanThenScore);
  if (VERB > 1) {
    std::cerr << "Before TDC there are " << scoreHolders.size() << " PSMs."
              << std::endl;
  }
  // First, sort the vector using the custom comparator

  // Use a map to track scan numbers and keep only the best score for each scan
  std::map<int, ScoreHolder*> bestScoreByScan;
  for (ScoreHolder* holder : scoreHolders) {
    int scan = holder->pPSM->scan;
    if (bestScoreByScan.find(scan) == bestScoreByScan.end() ||
        bestScoreByScan[scan]->score < holder->score) {
      bestScoreByScan[scan] = holder;
    }
  }

  // Clear the original scoreHolders and repopulate it with the best scores from
  // each scan
  scoreHolders.clear();
  for (const auto& pair : bestScoreByScan) {
    scoreHolders.push_back(pair.second);
  }
  if (VERB > 1) {
    std::cerr << "After TDC there are " << scoreHolders.size() << " PSMs."
              << std::endl;
  }
}

namespace {

// Builds the signature of generateCompositionSignature without a std::map
// per peptide: the residues (single characters) are counted in a table and
// the modifications in a short list. As std::string keys, the residues
// compared by their unsigned character, and the modifications, which start
// with '[', came after the residues below '[' and before those above it,
// and compared as strings among themselves; the signature lists them in
// that order.
class SignatureBuilder {
 public:
  void build(std::string_view peptide, std::string& signature) {
    for (size_t i = 0; i < peptide.size(); ++i) {
      if (peptide[i] == '[') {
        size_t j = i + 1;
        while (j < peptide.size() && peptide[j] != ']') {
          ++j;
        }
        // the modification, up to the end if it is not closed
        addModification(peptide.substr(i, j - i + 1));
        i = j;
      } else {
        const unsigned char residue = static_cast<unsigned char>(peptide[i]);
        if (counts_[residue]++ == 0) {
          residues_.push_back(residue);
        }
      }
    }
    std::sort(residues_.begin(), residues_.end());
    std::sort(modifications_.begin(), modifications_.end());
    signature.clear();
    size_t r = 0;
    for (; r < residues_.size() && residues_[r] < '['; ++r) {
      appendResidue(residues_[r], signature);
    }
    for (const auto& modification : modifications_) {
      signature += modification.first;
      signature += std::to_string(modification.second);
    }
    for (; r < residues_.size(); ++r) {
      appendResidue(residues_[r], signature);
    }
    for (unsigned char residue : residues_) {
      counts_[residue] = 0;
    }
    residues_.clear();
    modifications_.clear();
  }

 private:
  void addModification(std::string_view modification) {
    for (auto& counted : modifications_) {
      if (counted.first == modification) {
        ++counted.second;
        return;
      }
    }
    modifications_.emplace_back(modification, 1);
  }
  void appendResidue(unsigned char residue, std::string& signature) const {
    signature += static_cast<char>(residue);
    signature += std::to_string(counts_[residue]);
  }

  int counts_[256] = {};
  std::vector<unsigned char> residues_;
  std::vector<std::pair<std::string_view, int> > modifications_;
};

}  // namespace

std::string CompositionSorter::generateCompositionSignature(
    std::string_view peptide) {
  SignatureBuilder builder;
  std::string signature;
  builder.build(peptide, signature);
  return signature;
}

int CompositionSorter::addPSMs(const Scores& scores, bool useTDC) {
  std::vector<ScoreHolder*> scoreHolders;
  for (const auto& scr : scores) {
    scoreHolders.push_back(const_cast<ScoreHolder*>(&scr));
  }
  if (useTDC) {
    targetDecoyCompetition(scoreHolders);
  }
  SignatureBuilder builder;
  std::string signature;
  compositionPsms_.reserve(compositionPsms_.size() + scoreHolders.size());
  for (ScoreHolder* pScr : scoreHolders) {
    std::string_view peptide = pScr->getPSM()->getPeptideSequenceView();
    builder.build(peptide, signature);
    // a new signature gets the next number, so the compositions are
    // numbered in the order in which they first appear
    auto composition = compositionIndex_.find(signature);
    if (composition == compositionIndex_.end()) {
      composition =
          compositionIndex_
              .emplace(signature,
                       static_cast<unsigned int>(compositionIndex_.size()))
              .first;
    }
    compositionPsms_.push_back(CompositionPsm{
        composition->second, peptide, compositionPsms_.size(), pScr});
  }
  return 0;
}

int CompositionSorter::sortScorePerPeptide() {
  const size_t numCompositions = compositionIndex_.size();
  // Order the PSMs by the number of their composition, which is the order
  // in which the compositions first appeared, and not by the iteration
  // order of the hash map: that order depends on the hash function and on
  // the bucket layout of the standard library, so the same input and seed
  // gave other peptides on another platform. A counting sort keeps the
  // PSMs of a composition in the order in which they were added ...
  const size_t n = compositionPsms_.size();
  std::vector<size_t> compositionStart(numCompositions + 1, 0);
  for (const CompositionPsm& cp : compositionPsms_) {
    ++compositionStart[cp.composition + 1];
  }
  for (size_t r = 0; r < numCompositions; ++r) {
    compositionStart[r + 1] += compositionStart[r];
  }
  {
    // the new position of each PSM, applied in place, cycle by cycle
    std::vector<size_t> destination(n);
    std::vector<size_t> next(compositionStart.begin(),
                             compositionStart.end() - 1);
    for (size_t i = 0; i < n; ++i) {
      destination[i] = next[compositionPsms_[i].composition]++;
    }
    for (size_t i = 0; i < n; ++i) {
      while (destination[i] != i) {
        const size_t d = destination[i];
        std::swap(compositionPsms_[i], compositionPsms_[d]);
        std::swap(destination[i], destination[d]);
      }
    }
  }
  // ... and then by peptide within each composition, as the map of the
  // peptides of each composition iterated, with the PSMs of a peptide in
  // the order added
  auto byPeptide = [](const CompositionPsm& a, const CompositionPsm& b) {
    if (a.peptide != b.peptide) {
      return a.peptide < b.peptide;
    }
    return a.order < b.order;
  };
  // These sorts, and those of the PSMs of each peptide below, are small
  // apart from a few peptides with thousands of PSMs, so they run on the
  // calling thread: a team of threads for each of those costs more than it
  // saves. The comparators are total orders (for scores that are not NaN),
  // so the result does not depend on the sort algorithm.
  for (size_t r = 0; r < numCompositions; ++r) {
    std::sort(compositionPsms_.begin() + compositionStart[r],
              compositionPsms_.begin() + compositionStart[r + 1],
              byPeptide);
  }

  // Comparator for sorting ScoreHolder references
  auto compareScoreHolder = [](const CompositionPsm& lhs,
                               const CompositionPsm& rhs) {
    return *lhs.psm > *rhs.psm;  // descending score, ties in a total order
  };
  // sort the PSMs of each peptide
  for (size_t first = 0; first < n;) {
    size_t last = first + 1;
    while (last < n &&
           compositionPsms_[last].composition ==
               compositionPsms_[first].composition &&
           compositionPsms_[last].peptide == compositionPsms_[first].peptide) {
      ++last;
    }
    std::sort(compositionPsms_.begin() + first,
              compositionPsms_.begin() + last, compareScoreHolder);
    first = last;
  }
  return 0;
}

void incVector(std::vector<size_t>& sizes, size_t s) {
  // Ensure the vector is large enough to hold the element at index 's'
  if (s >= sizes.size()) {
    // Resize the vector to accommodate the new index, initializing new elements
    // to 0
    sizes.resize(s + 1, 0);
  }
  // Increment the element at index 's'
  sizes[s]++;
}

int CompositionSorter::inCompositionCompetition(Scores& bestScoreHolders,
                                                unsigned int decoysPerTarget) {
  // Sort the PSMs for each peptide in descending order of score
  sortScorePerPeptide();

  if (VERB > 1) {
    std::cerr << "Composition Matching starting with "
              << compositionIndex_.size() << " compositions." << std::endl;
  }

  std::vector<size_t> compSizeStat;
  std::vector<size_t> compTargetSizeStat;
  std::vector<size_t> compNumTargetSizeStat;
  const size_t numPsms = compositionPsms_.size();
  // the PSMs of a composition, and of each of its peptides, are adjacent
  for (size_t compFirst = 0, compLast = 0; compFirst < numPsms;
       compFirst = compLast) {
    const unsigned int composition = compositionPsms_[compFirst].composition;
    // the first PSM of each peptide of the composition
    std::vector<size_t> peptideFirsts;
    for (compLast = compFirst; compLast < numPsms &&
                               compositionPsms_[compLast].composition ==
                                   composition;
         ++compLast) {
      if (compLast == compFirst || compositionPsms_[compLast].peptide !=
                                       compositionPsms_[compLast - 1].peptide) {
        peptideFirsts.push_back(compLast);
      }
    }
    std::vector<std::vector<ScoreHolder*>> compositionGroups;
    std::vector<ScoreHolder*> targets;
    std::vector<ScoreHolder*> decoys;
    // register length statistics for printouts
    size_t numberOfPeptidesInComposition = peptideFirsts.size();
    incVector(compSizeStat, numberOfPeptidesInComposition);
    // Add the target peptides
    size_t numTargetPeptides = 0;
    for (size_t peptideFirst : peptideFirsts) {
      // This is the best scoring PSM for the peptide
      ScoreHolder* firstScoreHolder = compositionPsms_[peptideFirst].psm;
      if (firstScoreHolder->isTarget()) {
        targets.push_back(firstScoreHolder);
        numTargetPeptides++;
        incVector(compNumTargetSizeStat, numberOfPeptidesInComposition);
      } else {
        decoys.push_back(firstScoreHolder);
      }
    }
    if (numTargetPeptides > 0)
      incVector(compTargetSizeStat, numberOfPeptidesInComposition);
    // Format tuples of targets and decoys. Now we are not checking if there are
    // enough decoys for each target.
    // TODO: Check if there are enough decoys for each target
    std::reverse(decoys.begin(), decoys.end());
    for (ScoreHolder* target : targets) {
      std::vector<ScoreHolder*> group;
      group.push_back(target);  // Add the target ScoreHolder to the group

      // Add decoysPerTarget number of decoys to the group
      for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty(); ++i) {
        group.push_back(decoys.back());  // Add the last decoy (the one first
                                         // added to the vector)
        decoys.pop_back();  // Remove the added decoy from the decoys list
      }
      compositionGroups.push_back(
          group);  // Add the newly formed group to compositionGroups
    }
    // Take care of the remaining decoys
    for (ScoreHolder* decoy : decoys) {
      std::vector<ScoreHolder*> group;
      group.push_back(decoy);  // Add the target ScoreHolder to the group

      // Add decoysPerTarget number of decoys to the group
      for (unsigned int i = 1; i < decoysPerTarget && !decoys.empty(); ++i) {
        group.push_back(decoys.back());  // Add the last decoy (the one first
                                         // added to the vector)
        decoys.pop_back();  // Remove the added decoy from the decoys list
      }
      compositionGroups.push_back(
          group);  // Add the newly formed group to compositionGroups
    }
    // Sort each group in compositionGroups based on the ScoreHolder's score
    // and add it to bestScoreHolders
    for (auto& group : compositionGroups) {
      if (group.empty())
        continue;
      std::sort(group.begin(), group.end(),
                [](const ScoreHolder* a, const ScoreHolder* b) -> bool {
                  return *a > *b;  // descending score, ties in a total order
                });
      bestScoreHolders.addScoreHolder(
          *(group.front()));  // Add the highest-scoring ScoreHolder
    }
  }
  if (VERB > 1) {
    std::cerr
        << "Composition Group Statistics: Sizes of composition groups, sizes "
           "with at least 1 target peptide, and avg. num targets per group."
        << std::endl;
    std::cerr << "Size\tTotal\tIsTarget\tAvgTarget" << std::endl;
    for (size_t ix = compSizeStat.size(); --ix;) {
      std::cerr << ix << '\t' << compSizeStat[ix] << '\t';
      if (ix < compTargetSizeStat.size()) {
        std::cerr << compTargetSizeStat[ix] << '\t'
                  << compTargetSizeStat[ix] / ix << std::endl;
      } else {
        std::cerr << 0 << '\t' << 0 << std::endl;
      }
    }
  }

  // counted here: Scores::size() comes from the totals, which the caller
  // only updates afterwards with recalculateSizes()
  int numPeptides = 0, numTargets = 0;
  for (const ScoreHolder& pSH : bestScoreHolders) {
    numPeptides++;
    if (pSH.isTarget())
      numTargets++;
  }
  if (VERB > 1) {
    std::cerr << "Composition Matching ends with " << numPeptides
              << " peptides, whereof " << numTargets << " is target peptides."
              << std::endl;
  }
  return 0;
}

int CompositionSorter::psmAndPeptide(const Scores& scores,
                                     Scores& winnerPeptides,
                                     unsigned int decoysPerTarget) {
  //    psmLevelCompetition(); // Should be handled by the reader?

  // Populate the structure with the winning PSMs
  addPSMs(scores);

  // Split out tuples of peptides of identical composition, and select the most
  // high scoring peptide in each tuple
  inCompositionCompetition(winnerPeptides, decoysPerTarget);
  winnerPeptides.recalculateSizes();
  return 0;
}

void CompositionSorter::psmsOnly(const Scores& scores, Scores& winnerPeptides) {
  // The best ScoreHolder of each peptide, in the order in which the
  // peptides first appear in scores. The hash map only looks up the
  // position of a peptide: the winners are not read in its iteration order,
  // which depends on std::hash and the bucket layout of the standard
  // library and so differed between libstdc++, libc++ and MSVC. The keys
  // view the peptides of the PSMs, which outlive the map.
  std::unordered_map<std::string_view, std::size_t> peptideIndex;
  std::vector<const ScoreHolder*> bestScoreHolders;

  // Iterate over each ScoreHolder, with one lookup per PSM
  for (const ScoreHolder& sh : scores) {
    auto inserted = peptideIndex.try_emplace(
        sh.getPSM()->getPeptideSequenceView(), bestScoreHolders.size());
    if (inserted.second) {
      bestScoreHolders.push_back(&sh);
    } else if (sh.score > bestScoreHolders[inserted.first->second]->score) {
      // the current ScoreHolder has a higher score for the peptide
      bestScoreHolders[inserted.first->second] = &sh;
    }
  }

  for (const ScoreHolder* best : bestScoreHolders) {
    winnerPeptides.addScoreHolder(*best);
  }
  winnerPeptides.recalculateSizes();

  if (VERB > 1) {
    std::cerr << "Selected the best scoring PSM for each of the "
              << winnerPeptides.size() << " peptides from a dataset of "
              << scores.size() << " PSMs." << std::endl;
  }
}

void CompositionSorter::retainRepresentatives(const Scores& psms,
                                              Scores& winnerPeptides,
                                              double /* selectionFDR */,
                                              unsigned int decoysPerTarget,
                                              bool useCompositionMatch) {
  if (useCompositionMatch) {
    if (VERB > 1) {
      std::cerr << "Starting reset: psmAndPeptide" << std::endl;
    }
    CompositionSorter sorter;
    sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);
  } else {
    if (VERB > 1) {
      std::cerr << "Starting reset: psmsOnly" << std::endl;
    }
    psmsOnly(psms, winnerPeptides);
  }
}