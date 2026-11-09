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
#ifndef PSM_STORE_H_
#define PSM_STORE_H_

#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "PSMDescription.h"

/*
 * PsmStore owns PSMDescription objects in fixed-size blocks, together with
 * the characters of their ids and peptides and their protein lists. All
 * addresses stay valid while more PSMs are added, and there is no heap
 * allocation per PSM. Protein names are stored once each.
 */
class PsmStore {
 public:
  PsmStore() = default;
  PsmStore(const PsmStore&) = delete;
  PsmStore& operator=(const PsmStore&) = delete;

  // moves psm into the store and returns its address in the store
  PSMDescription* add(PSMDescription&& psm) {
    PSMDescription* slot;
    if (!freeSlots_.empty()) {
      slot = freeSlots_.back();
      freeSlots_.pop_back();
    } else {
      if (numUsed_ == blocks_.size() * kBlockSize) {
        blocks_.emplace_back(new PSMDescription[kBlockSize]);
      }
      slot = &blocks_.back()[numUsed_ % kBlockSize];
      ++numUsed_;
    }
    *slot = std::move(psm);
    return slot;
  }

  // empties the slot of a PSM so that a later add can reuse it
  void release(PSMDescription* psm) {
    *psm = PSMDescription();
    freeSlots_.push_back(psm);
  }

  // copies the characters into the store and returns a view of the copy
  std::string_view storeString(std::string_view s) {
    if (s.size() > kCharBlockSize / 4) {  // too long to share a block
      longStrings_.emplace_back(new char[s.size()]);
      std::memcpy(longStrings_.back().get(), s.data(), s.size());
      return std::string_view(longStrings_.back().get(), s.size());
    }
    if (charBlocks_.empty() || charUsed_ + s.size() > kCharBlockSize) {
      charBlocks_.emplace_back(new char[kCharBlockSize]);
      charUsed_ = 0;
    }
    char* dest = charBlocks_.back().get() + charUsed_;
    std::memcpy(dest, s.data(), s.size());
    charUsed_ += s.size();
    return std::string_view(dest, s.size());
  }

  // stores a list of protein names, each name only once in the store, and
  // returns the first element of the stored list of n names
  const std::string* const* storeProteins(
      const std::vector<std::string_view>& names) {
    const size_t n = names.size();
    if (n == 0) return nullptr;
    const std::string** list;
    if (n > kListBlockSize / 4) {  // too long to share a block
      longLists_.emplace_back(new const std::string*[n]);
      list = longLists_.back().get();
    } else {
      if (listBlocks_.empty() || listUsed_ + n > kListBlockSize) {
        listBlocks_.emplace_back(new const std::string*[kListBlockSize]);
        listUsed_ = 0;
      }
      list = listBlocks_.back().get() + listUsed_;
      listUsed_ += n;
    }
    for (size_t i = 0; i < n; ++i) {
      auto found = proteinIndex_.find(names[i]);
      if (found == proteinIndex_.end()) {
        proteinNames_.emplace_back(names[i]);  // deque elements never move
        found = proteinIndex_
                    .emplace(std::string_view(proteinNames_.back()),
                             &proteinNames_.back())
                    .first;
      }
      list[i] = found->second;
    }
    return list;
  }

  void clear() {
    blocks_.clear();
    freeSlots_.clear();
    numUsed_ = 0;
    charBlocks_.clear();
    longStrings_.clear();
    charUsed_ = 0;
    listBlocks_.clear();
    longLists_.clear();
    listUsed_ = 0;
    proteinIndex_.clear();
    proteinNames_.clear();
  }

 private:
  static const size_t kBlockSize = 1u << 16;
  std::vector<std::unique_ptr<PSMDescription[]> > blocks_;
  std::vector<PSMDescription*> freeSlots_;
  size_t numUsed_ = 0;

  // strings and lists are packed into the last block; long ones get their own
  static const size_t kCharBlockSize = 1u << 20;
  std::vector<std::unique_ptr<char[]> > charBlocks_, longStrings_;
  size_t charUsed_ = 0;
  static const size_t kListBlockSize = 1u << 16;
  std::vector<std::unique_ptr<const std::string*[]> > listBlocks_, longLists_;
  size_t listUsed_ = 0;
  std::deque<std::string> proteinNames_;
  std::unordered_map<std::string_view, const std::string*> proteinIndex_;
};

#endif /* PSM_STORE_H_ */
