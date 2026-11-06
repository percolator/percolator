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
 * FASTA reading, in-silico digestion and a peptide-to-protein index for the
 * detection of fragment and duplicate proteins (PickedProteinCaller).
 */
#ifndef PICKED_PROTEIN_PROTEIN_DIGEST_H_
#define PICKED_PROTEIN_PROTEIN_DIGEST_H_

#include <cstddef>
#include <cstdint>
#include <istream>
#include <string>
#include <vector>

class Enzyme;

namespace PickedProtein {

/** Which peptide ends must be cleavage sites of the enzyme. */
enum class DigestType {
  Full,        ///< both ends (a protein terminus counts as a cleavage site)
  Partial,     ///< at least one end
  NonSpecific  ///< neither end
};

struct DigestConstraint {
  DigestType digest;
  int minLength;
  int maxLength;
  /** Cleavage sites allowed inside a peptide; NonSpecific digests allow 500. */
  int maxMissedCleavages;
};

/** A protein of a FASTA file. */
struct FastaProtein {
  std::string id;        ///< first word of the header line
  std::string sequence;  ///< letters of the sequence, in upper case
};

/**
 * Appends the proteins of a FASTA stream to proteins. Lines before the first
 * line that starts with '>' are skipped. After that, every '>' starts a
 * header, which extends to the end of its line; the sequence is everything up
 * to the next '>' or the end of the stream, of which only the letters A-Z and
 * a-z are kept, in upper case. There is no length limit.
 * Returns false if reading failed (badbit).
 */
bool readFasta(std::istream& in, std::vector<FastaProtein>& proteins);

/** A peptide as a 0-based start and a length in its protein's sequence. */
struct PeptideSpan {
  uint32_t start;
  uint32_t length;
};

/**
 * Appends the peptides of sequence that satisfy constraint to spans, in an
 * unspecified order. The cleavage sites are those of Enzyme::isEnzymatic;
 * the residue after the C-terminal one is given as '\0', so for it the rule
 * decides on the C-terminal residue alone. Missed cleavages are the
 * cleavage sites strictly inside the peptide. Partial digests reproduce the
 * PercolatorCrux digestion this replaces, which left out the peptides that
 * start at a cleavage site and end at the last non-cleavage position, when
 * that position is not the C terminus.
 */
void digestProtein(const std::string& sequence, const Enzyme& enzyme,
                   const DigestConstraint& constraint,
                   std::vector<PeptideSpan>& spans);

/**
 * Peptide-to-protein index of a digest of some proteins of a database. Each
 * peptide that starts at the protein N terminus with a methionine is also
 * indexed without that methionine if it then still has minLength residues.
 * A peptide that occurs several times in a protein is indexed as often.
 * The occurrences are kept in flat arrays grouped by sequence, about 8 bytes
 * per occurrence; building them needs 24 more bytes per occurrence.
 */
class PeptideIndex {
 public:
  /** Indexes proteins[i] for each i of proteinIds, which must be ascending;
      an empty proteinIds stands for all proteins. */
  PeptideIndex(const std::vector<FastaProtein>& proteins,
               const std::vector<size_t>& proteinIds, const Enzyme& enzyme,
               const DigestConstraint& constraint);

  /** Peptides of an indexed protein, without the methionine-cleaved ones. */
  size_t numPeptides(size_t proteinIdx) const;

  /**
   * Sets result to the proteins that contain every indexed peptide of
   * proteinIdx, the protein itself included: the intersection of the
   * ascending protein lists of its peptides, counting a protein as often as
   * it has the peptide. Stops as soon as fewer than two proteins remain, so
   * a result with fewer than two proteins may be incomplete. Empty if the
   * protein has no peptides.
   */
  void proteinsContainingAllPeptides(size_t proteinIdx,
                                     std::vector<size_t>& result) const;

 private:
  size_t position(size_t proteinIdx) const;

  std::vector<size_t> proteinIds_;  // empty: all proteins, position = index
  std::vector<size_t> numPeptides_;
  // occurrences of the protein at position k, in digestion order:
  // [occurrenceBegin_[k], occurrenceBegin_[k + 1])
  std::vector<size_t> occurrenceBegin_;
  std::vector<uint32_t> groupOfOccurrence_;
  // proteins of the occurrences of group g, ascending:
  // groupProteins_[groupBegin_[g] .. groupBegin_[g + 1])
  std::vector<uint32_t> groupBegin_;
  std::vector<uint32_t> groupProteins_;
};

}  // namespace PickedProtein

#endif /* PICKED_PROTEIN_PROTEIN_DIGEST_H_ */
