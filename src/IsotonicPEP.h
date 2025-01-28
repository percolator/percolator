/*******************************************************************************
 Copyright 2024 Lukas Käll <lukas.kall@scilifelab.se>

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

#ifndef ISOTONICPEP_H_
#define ISOTONICPEP_H_

#include <iostream>
#include <vector>
#include <algorithm>  // for std::sort, std::transform, std::adjacent_difference
#include <numeric>    // for std::iota
#include <stdexcept>
#include <cmath>      // for std::abs


class IsotonicPEP {

protected:
    struct Block {
        double sum;
        int count;
        double avg;
    };

public:
    IsotonicPEP() = default;

    // Method to perform (non-decreasing) isotonic regression

    std::vector<double> isotonic_regression_pava(const std::vector<double>& y) {
        const int n = static_cast<int>(y.size());
        if (n == 0) {
            return {};
        }

        // We'll store our result here
        std::vector<double> result(n, 0.0);

        // Stack of "blocks", each block represents a contiguous segment
        // that has been merged so far
        std::vector<Block> stack;
        stack.reserve(n);

        // 1. Go through each data point (left to right)
        for (int i = 0; i < n; ++i) {
            // Create a block for the single point y[i]
            Block newBlock { y[i], 1, y[i] };
            stack.push_back(newBlock);

            // 2. While there's a violation of the non-decreasing property
            //    (the previous block's average > current block's average),
            //    merge the top two blocks
            while (stack.size() > 1) {
                // top block
                Block &top     = stack[stack.size() - 1];
                // second-to-top block
                Block &secTop  = stack[stack.size() - 2];

                if (secTop.avg > top.avg) {
                    // Merge them
                    double mergedSum   = secTop.sum + top.sum;
                    int    mergedCount = secTop.count + top.count;
                    double mergedAvg   = mergedSum / mergedCount;

                    // Pop the two blocks off
                    stack.pop_back();
                    stack.pop_back();

                    // Push the merged block
                    Block mergedBlock { mergedSum, mergedCount, mergedAvg };
                    stack.push_back(mergedBlock);
                } else {
                    // No violation: non-decreasing order is satisfied
                    break;
                }
            }
        }

        // 3. Now expand the merged blocks into the final piecewise-constant solution
        int index = 0;
        for (const auto& block : stack) {
            for (int c = 0; c < block.count; ++c) {
                result[index++] = block.avg;
            }
        }
        // Clip the values to be within the bounds [0, 1]
        for (auto& val : result) {
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
        }
 
        return result;
    }


    std::vector<double> q_to_pep(const std::vector<double>& q_values) {
        qs = q_values;

        // Generate the qn values as q-values multiplied by their respective indices
        std::vector<double> qn(q_values.size());
        std::vector<int> indices(q_values.size());
        std::iota(indices.begin(), indices.end(), 1);  // Generate indices starting from 1

        for (size_t i = 0; i < q_values.size(); ++i) {
            qn[i] = q_values[i] * indices[i];
            assert((i == q_values.size()-1) || (q_values[i] <=  q_values[i+1]) );
        }

        // Calculate differences between consecutive qn values

        std::vector<double> raw_pep(qn.size());
        raw_pep[0] = qn[0];
        for (size_t i = 1; i < qn.size(); ++i) {
            raw_pep[i] = qn[i] - qn[i-1];
        }
        
        // Perform isotonic regression on the differences
        std::vector<double> pep_iso = isotonic_regression_pava(raw_pep);
//        for (size_t i = 0; i < qn.size()-1; ++i) {
//            cerr << i << " " << q_values[i] << " " << qn[i] << " " << raw_pep[i] << " " << pep_iso[i] << endl;
//        }
        return pep_iso;
    }

    double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
        double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
        return interp_pep;
    }

    // Getter method to interpolate and retrieve PEP for a given q-value
    double get_pep(double q_value) const {
        if (qs.empty() || pep_iso.empty()) {
            throw std::runtime_error("PEP calculation has not been performed yet.");
        }

        auto it = std::lower_bound(qs.begin(), qs.end(), q_value);

        if (it == qs.begin()) {
            // If the q_value is less than or equal to the smallest q, return the first PEP
            return pep_iso.front();
        } else if (it == qs.end()) {
            // If the q_value is greater than the largest q, return the last PEP
            return pep_iso.back();
        } else {
            // Otherwise, interpolate between the two closest points
            size_t idx = std::distance(qs.begin(), it);
            double q1 = qs[idx - 1];
            double q2 = qs[idx];
            double pep1 = pep_iso[idx - 1];
            double pep2 = pep_iso[idx];

            // Linear interpolation
            double interp_pep = interpolate(q_value, q1, q2, pep1, pep2);
            return interp_pep;
        }
    }

protected:
    std::vector<double> qs;
    std::vector<double> pep_iso;

};
#endif /* ISOTONICPEP_H_ */
