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
#ifndef SCAN_ID_TABLE_H_
#define SCAN_ID_TABLE_H_

#include <algorithm>
#include <cstddef>
#include <vector>

/*
 * Tells whether a scan id, i.e. a scan number and an experimental mass,
 * occurs both in a PSM with a decoy label and in a PSM with another label,
 * which means that the targets and decoys come from separate searches.
 *
 * Scan ids are equal when both the scan numbers and the masses compare
 * equal, as in a hash table keyed on the pair: 0.0 and -0.0 are the same
 * mass, and a NaN mass is never equal to another one.
 *
 * The scan ids are collected in a vector, which costs 16 bytes per PSM
 * instead of a hash table node, and sorted. The vector is sorted and checked
 * whenever its size doubles, by sorting the new part and merging it into the
 * sorted part. Separate searches therefore usually show up early, after which
 * the vector is released and further scan ids are ignored.
 */
class SeparateSearchDetector {
 public:
  // returns true once a scan id with both labels has been found
  bool add(int scan, double mass, bool isDecoy) {
    if (found_) return true;
    if (mass != mass) return false;  // a NaN mass matches no other scan id
    keys_.push_back(Key{scan, isDecoy, mass});
    if (keys_.size() >= nextCheck_) {
      check();
      nextCheck_ = 2 * keys_.size();
    }
    return found_;
  }

  // checks the scan ids added since the last check and releases the memory
  bool foundSeparateSearches() {
    if (!found_ && sorted_ < keys_.size()) check();
    std::vector<Key>().swap(keys_);
    sorted_ = 0;
    return found_;
  }

 private:
  struct Key {
    int scan;
    bool isDecoy;
    double mass;
  };
  // Orders by scan id, then by label. A total order on scan ids without
  // NaN masses: 0.0 and -0.0 are equivalent, as they are equal.
  static bool less(const Key& a, const Key& b) {
    if (a.scan != b.scan) return a.scan < b.scan;
    if (a.mass != b.mass) return a.mass < b.mass;
    return a.isDecoy < b.isDecoy;
  }

  void check() {
    std::vector<Key>::iterator mid = keys_.begin() + static_cast<std::ptrdiff_t>(sorted_);
    std::sort(mid, keys_.end(), less);
    std::inplace_merge(keys_.begin(), mid, keys_.end(), less);
    sorted_ = keys_.size();
    // within a run of equal scan ids, a decoy follows a non-decoy
    for (size_t i = 1; i < keys_.size(); ++i) {
      if (keys_[i].isDecoy != keys_[i - 1].isDecoy &&
          keys_[i].scan == keys_[i - 1].scan &&
          keys_[i].mass == keys_[i - 1].mass) {
        found_ = true;
        std::vector<Key>().swap(keys_);
        sorted_ = 0;
        return;
      }
    }
  }

  std::vector<Key> keys_;
  size_t sorted_ = 0;  // keys_[0, sorted_) is sorted
  size_t nextCheck_ = 1u << 16;
  bool found_ = false;
};

#endif /* SCAN_ID_TABLE_H_ */
