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

#include "ProteinDigest.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <string_view>

#include "Enzyme.h"
#include "MyException.h"

namespace PickedProtein {

namespace {

bool isBlank(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' ||
         c == '\r';
}

std::string firstWord(const std::string& header) {
  size_t begin = 0;
  while (begin < header.size() && isBlank(header[begin])) ++begin;
  size_t end = begin;
  while (end < header.size() && !isBlank(header[end])) ++end;
  return header.substr(begin, end - begin);
}

void appendResidues(const std::string& text, std::string& sequence) {
  for (char c : text) {
    if (c >= 'A' && c <= 'Z') {
      sequence.push_back(c);
    } else if (c >= 'a' && c <= 'z') {
      sequence.push_back(static_cast<char>(c - 'a' + 'A'));
    }
  }
}

}  // namespace

bool readFasta(std::istream& in, std::vector<FastaProtein>& proteins) {
  std::string header;
  bool inProtein = false;
  while (std::getline(in, header)) {
    if (!header.empty() && header[0] == '>') {
      header.erase(0, 1);
      inProtein = true;
      break;
    }
  }
  std::string text;
  while (inProtein) {
    FastaProtein protein;
    protein.id = firstWord(header);
    text.clear();  // getline leaves it unchanged on a failed stream
    std::getline(in, text, '>');
    protein.sequence.reserve(text.size());  // no regrowth while appending
    appendResidues(text, protein.sequence);
    proteins.push_back(std::move(protein));
    // a '>' was found unless the stream ended
    inProtein = in.good();
    if (inProtein) std::getline(in, header);
  }
  return !in.bad();
}

void digestProtein(const std::string& sequence, const Enzyme& enzyme,
                   const DigestConstraint& constraint,
                   std::vector<PeptideSpan>& spans) {
  const size_t length = sequence.size();
  if (length == 0) return;
  // isCut[p]: position p (between residues p-1 and p) is a cleavage site;
  // the protein termini count as cleavage sites.
  std::vector<char> isCut(length + 1);
  isCut[0] = 1;
  size_t lastNonCut = 0;  // 0: none, as position 0 is always a cut
  for (size_t i = 0; i < length; ++i) {
    const char next = i + 1 < length ? sequence[i + 1] : '\0';
    if (enzyme.isEnzymatic(sequence[i], next)) {
      isCut[i + 1] = 1;
    } else {
      lastNonCut = i + 1;
    }
  }
  isCut[length] = 1;

  const DigestType digest = constraint.digest;
  const long maxMissed = digest == DigestType::NonSpecific
                             ? 500
                             : constraint.maxMissedCleavages;
  const long minLength = constraint.minLength;
  const long maxLength = constraint.maxLength;
  if (maxLength < 1) return;
  for (size_t start = 0; start < length; ++start) {
    if (digest == DigestType::Full && !isCut[start]) continue;
    const size_t endMax =
        std::min(length, start + static_cast<size_t>(maxLength));
    long missed = 0;  // cleavage sites strictly between start and end
    for (size_t end = start + 1; end <= endMax; ++end) {
      if (missed > maxMissed) break;
      if (static_cast<long>(end - start) >= minLength) {
        bool keep = true;
        if (digest == DigestType::Full) {
          keep = isCut[end];
        } else if (digest == DigestType::Partial) {
          keep = isCut[end] || (isCut[start] && end != lastNonCut);
        }
        if (keep) {
          spans.push_back({static_cast<uint32_t>(start),
                           static_cast<uint32_t>(end - start)});
        }
      }
      if (isCut[end]) ++missed;
    }
  }
}

namespace {

struct Occurrence {
  uint64_t hash;
  uint32_t protein;
  uint32_t start;
  uint32_t length;
  uint32_t index;  // position in digestion order
};

}  // namespace

PeptideIndex::PeptideIndex(const std::vector<FastaProtein>& proteins,
                           const std::vector<size_t>& proteinIds,
                           const Enzyme& enzyme,
                           const DigestConstraint& constraint)
    : proteinIds_(proteinIds) {
  const size_t numIndexed =
      proteinIds_.empty() ? proteins.size() : proteinIds_.size();
  if (proteins.size() > std::numeric_limits<uint32_t>::max()) {
    throw MyException("ERROR: too many proteins in the FASTA database.");
  }
  numPeptides_.resize(numIndexed);
  occurrenceBegin_.resize(numIndexed + 1);

  std::vector<PeptideSpan> spans;
  // count the occurrences first, so that their array is allocated once
  size_t numOccurrences = 0;
  for (size_t k = 0; k < numIndexed; ++k) {
    const size_t protein = proteinIds_.empty() ? k : proteinIds_[k];
    const std::string& sequence = proteins[protein].sequence;
    spans.clear();
    digestProtein(sequence, enzyme, constraint, spans);
    numOccurrences += spans.size();
    for (const PeptideSpan& span : spans) {
      if (span.start == 0 && sequence[0] == 'M' &&
          static_cast<long>(span.length) - 1 >= constraint.minLength) {
        ++numOccurrences;
      }
    }
  }
  if (numOccurrences > std::numeric_limits<uint32_t>::max()) {
    throw MyException(
        "ERROR: too many peptides in the digest of the FASTA database.");
  }
  std::vector<Occurrence> occurrences;
  occurrences.reserve(numOccurrences);
  std::hash<std::string_view> hasher;
  for (size_t k = 0; k < numIndexed; ++k) {
    const size_t protein = proteinIds_.empty() ? k : proteinIds_[k];
    const std::string& sequence = proteins[protein].sequence;
    const std::string_view view(sequence);
    occurrenceBegin_[k] = occurrences.size();
    spans.clear();
    digestProtein(sequence, enzyme, constraint, spans);
    numPeptides_[k] = spans.size();
    for (const PeptideSpan& span : spans) {
      occurrences.push_back({hasher(view.substr(span.start, span.length)),
                             static_cast<uint32_t>(protein), span.start,
                             span.length, 0});
      if (span.start == 0 && sequence[0] == 'M' &&
          static_cast<long>(span.length) - 1 >= constraint.minLength) {
        occurrences.push_back({hasher(view.substr(1, span.length - 1)),
                               static_cast<uint32_t>(protein), 1,
                               span.length - 1, 0});
      }
    }
  }
  occurrenceBegin_[numIndexed] = numOccurrences;
  for (size_t i = 0; i < numOccurrences; ++i) {
    occurrences[i].index = static_cast<uint32_t>(i);
  }

  auto sequenceOf = [&proteins](const Occurrence& o) {
    return std::string_view(proteins[o.protein].sequence)
        .substr(o.start, o.length);
  };
  // group equal sequences, each group's proteins ascending
  auto bySequence = [&sequenceOf](const Occurrence& a, const Occurrence& b) {
    if (a.hash != b.hash) return a.hash < b.hash;
    const int cmp = sequenceOf(a).compare(sequenceOf(b));
    if (cmp != 0) return cmp < 0;
    return a.protein < b.protein;
  };
#if defined(_GLIBCXX_PARALLEL) && defined(__GLIBCXX__)
  // the parallel mergesort of libstdc++ would need a copy of all occurrences
  std::sort(occurrences.begin(), occurrences.end(), bySequence,
            __gnu_parallel::sequential_tag());
#else
  std::sort(occurrences.begin(), occurrences.end(), bySequence);
#endif
  groupOfOccurrence_.resize(numOccurrences);
  groupProteins_.resize(numOccurrences);
  for (size_t i = 0; i < numOccurrences; ++i) {
    const Occurrence& o = occurrences[i];
    if (i == 0 || o.hash != occurrences[i - 1].hash ||
        sequenceOf(o) != sequenceOf(occurrences[i - 1])) {
      groupBegin_.push_back(static_cast<uint32_t>(i));
    }
    groupProteins_[i] = o.protein;
    groupOfOccurrence_[o.index] = static_cast<uint32_t>(groupBegin_.size() - 1);
  }
  groupBegin_.push_back(static_cast<uint32_t>(numOccurrences));
}

size_t PeptideIndex::position(size_t proteinIdx) const {
  if (proteinIds_.empty()) return proteinIdx;
  return static_cast<size_t>(
      std::lower_bound(proteinIds_.begin(), proteinIds_.end(), proteinIdx) -
      proteinIds_.begin());
}

size_t PeptideIndex::numPeptides(size_t proteinIdx) const {
  return numPeptides_[position(proteinIdx)];
}

void PeptideIndex::proteinsContainingAllPeptides(
    size_t proteinIdx, std::vector<size_t>& result) const {
  const size_t k = position(proteinIdx);
  result.clear();
  std::vector<size_t> intersection;
  for (size_t i = occurrenceBegin_[k]; i < occurrenceBegin_[k + 1]; ++i) {
    const uint32_t group = groupOfOccurrence_[i];
    const uint32_t* first = groupProteins_.data() + groupBegin_[group];
    const uint32_t* last = groupProteins_.data() + groupBegin_[group + 1];
    if (i == occurrenceBegin_[k]) {
      result.assign(first, last);
    } else {
      intersection.clear();
      std::set_intersection(result.begin(), result.end(), first, last,
                            std::back_inserter(intersection));
      result.swap(intersection);
    }
    if (result.size() < 2) break;
  }
}

}  // namespace PickedProtein
