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

class IsotonicRegression {
public:
    IsotonicRegression() = default;

    // Method to perform isotonic regression
    std::vector<double> fit_transform(const std::vector<double>& y) {
        if (y.empty()) {
            throw std::invalid_argument("Input vector 'y' is empty.");
        }

        std::vector<double> solution = y;
        std::vector<int> level_start(y.size());

        for (size_t i = 0; i < y.size(); ++i) {
            level_start[i] = i;
        }

        for (size_t i = 1; i < solution.size(); ++i) {
            if (solution[i] < solution[i - 1]) {
                int j = i;
                while (j > 0 && solution[level_start[j - 1]] > solution[i]) {
                    --j;
                }
                double sum = 0.0;
                int start = level_start[j];
                int count = i - start + 1;
                for (int k = start; k <= i; ++k) {
                    sum += solution[k];
                }
                double average = sum / count;
                for (int k = start; k <= i; ++k) {
                    solution[k] = average;
                }
                level_start[i] = start;
            }
        }

        // Clip the values to be within the bounds [0, 1]
        for (auto& val : solution) {
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
        }

        return solution;
    }

    std::vector<double> q_to_pep(const std::vector<double>& q_values) {
        qs = q_values;

        // Generate the qn values as q-values multiplied by their respective indices
        std::vector<double> qn(q_values.size());
        std::vector<int> indices(q_values.size());
        std::iota(indices.begin(), indices.end(), 1);  // Generate indices starting from 1

        for (size_t i = 0; i < q_values.size(); ++i) {
            qn[i] = q_values[i] * indices[i];
        }

        // Calculate differences between consecutive qn values
        std::vector<double> raw_pep(qn.size() - 1);
        std::adjacent_difference(qn.begin(), qn.end(), raw_pep.begin());
        raw_pep.erase(raw_pep.begin());  // Remove the first element (which is meaningless)

        // Append an additional element (1.0) to match the original array size
        raw_pep.push_back(1.0);

        // Perform isotonic regression on the differences
        std::vector<double> pep_iso = fit_transform(raw_pep);
        return pep_iso;
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
            double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
            return interp_pep;
        }
    }

protected:
    std::vector<double> qs;
    std::vector<double> pep_iso;

};
#endif /* ISOTONICPEP_H_ */
