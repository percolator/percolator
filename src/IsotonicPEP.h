/*******************************************************************************
 Copyright 2025 Lukas Käll <lukas.kall@scilifelab.se>

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
#include <limits>
#include <iostream>


/**
 * Base class for isotonic regression in real space (no logistic transformations).
 * Provides:
 *   - pavaNonDecreasing(...) => Stepwise constant
 *   - pavaNonDecreasingInterpolation(...) => Linear interpolation between block centers
 */
class IsotonicRegression
{
protected:
    //--------------------------------------------------------
    // A simple block structure for merging
    //--------------------------------------------------------
    struct Block {
        double sum   = 0.0;
        int    count = 0;
        double avg   = 0.0;
    };

    //--------------------------------------------------------
    // An "extended" block that also tracks the range of
    // indices covered, for interpolation usage.
    //--------------------------------------------------------
    struct ExtBlock {
        double sum      = 0.0;
        int    count    = 0;
        double avg      = 0.0;
        int    startIdx = 0;  // inclusive
        int    endIdx   = 0;  // inclusive
    };

public:
    IsotonicRegression() = default;
    virtual ~IsotonicRegression() = default;

    /**
     * Perform standard PAVA (pool-adjacent-violators) to enforce a non-decreasing sequence
     * on the input 'values'.  
     * Here we extended the method to also merge block if average is outside 
     * the range [min_value, max_value] 
     *
     * Each final "block" is repeated block.count times in the output.
     *
     * Typically used in 1D for y-values only.
     */
    virtual std::vector<double> pavaNonDecreasingRanged(
        const std::vector<double>& values,
        const double min_value = std::numeric_limits<double>::min(),
        const double max_value = std::numeric_limits<double>::max()
    ) const
    {
        const int n = static_cast<int>(values.size());
        if (n == 0) {
            return {};
        }

        // We'll store blocks on a stack
        struct MergeBlock {
            double sum;
            int    count;
            double avg;
        };

        std::vector<MergeBlock> stack;
        stack.reserve(n);

        // 1. Left to right
        for (int i = 0; i < n; ++i) {
            MergeBlock newBlock { values[i], 1, values[i] };
            stack.push_back(newBlock);

            // 2. Merge while there's a violation of non-decreasing
            while (stack.size() > 1) {
                auto &top    = stack.back();
                auto &secTop = stack[stack.size() - 2];
                double mergedSum   = secTop.sum + top.sum;
                int    mergedCount = secTop.count + top.count;
                double mergedAvg   = mergedSum / mergedCount;

                // if (( secTop.avg > top.avg) || ( mergedAvg < min_value ) || ( mergedAvg > max_value )) {
                if ( secTop.avg > top.avg) {
                    stack.pop_back();
                    stack.pop_back();
                    stack.push_back({ mergedSum, mergedCount, mergedAvg });
                } else {
                    break;
                }
            }
        }

        // 3. Expand final solution
        std::vector<double> result;
        result.reserve(n); // approximate; actual size = sum(counts)

        for (auto &b : stack) {
            double val = max( min(b.avg, max_value), min_value);
            for (int c = 0; c < b.count; ++c) {
                result.push_back(val);
            }
        }
        return result;
    }


    /**
     * "Interpolated" variant of non-decreasing PAVA.
     * Instead of returning stepwise-constant blocks,
     * it returns a piecewise-linear function across block centers.
     *
     * For this, you need:
     *   - x[i]: sorted "positions"
     *   - y[i]: values to be made non-decreasing
     *
     * We'll:
     *   (a) Group the data into extended blocks (each is 1 point by default),
     *       then apply standard PAVA merges in non-decreasing order. So the
     *       final stack is a set of blocks. Each block covers [startIdx..endIdx]
     *       with a final average "avg".
     *   (b) For each final block, compute the "center" in x (could be midpoint or
     *       average of x). We'll store that in ExtBlock.
     *   (c) For each point i in that block, linearly interpolate in [center(block), center(nextBlock)]
     *       from avg(block) to avg(nextBlock).
     *
     * This does *not* strictly minimize the usual isotonic objective, but
     * might give you smoother transitions.
     */
    virtual std::vector<double> pavaNonDecreasingInterpolationRanged(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const double min_y = std::numeric_limits<double>::min(),
        const double max_y = std::numeric_limits<double>::max()
    ) const
    {
        if (x.size() != y.size()) {
            throw std::invalid_argument("pavaNonDecreasingInterpolation: x.size() != y.size()");
        }
        const int n = static_cast<int>(y.size());
        if (n == 0) {
            return {};
        }

        // 1. Build an ExtBlock for each point initially
        std::vector<ExtBlock> blocks;
        blocks.reserve(n);
        for (int i = 0; i < n; ++i) {
            blocks.push_back(ExtBlock {
                y[i],       // sum
                1,          // count
                y[i],       // avg
                i,          // startIdx
                i           // endIdx
            });
        }

        // 2. Standard PAVA merges in non-decreasing order
        //    We'll keep them in "blocks" but get a new final stack.
        //    Here we extended the method to also merge block if average is outside 
        //    the range [min_value, max_value] 
        std::vector<ExtBlock> stack;
        stack.reserve(n);

        for (auto &b : blocks) {
            stack.push_back(b);

            // Merge while there's a violation or out of range
            while (stack.size() > 1) {
                auto &top    = stack.back();
                auto &secTop = stack[stack.size()-2];
                double mergedSum   = secTop.sum + top.sum;
                int    mergedCount = secTop.count + top.count;
                double mergedAvg   = mergedSum / mergedCount;
                // if (( secTop.avg > top.avg) || ( mergedAvg < min_y ) || ( mergedAvg > max_y )) {
                if ( secTop.avg > top.avg ) {
                    int    sIdx     = secTop.startIdx;
                    int    eIdx     = top.endIdx;
                    stack.pop_back();
                    stack.pop_back();
                    stack.push_back({
                        mergedSum, mergedCount, mergedAvg, sIdx, eIdx
                    });
                } else {
                    break;
                }
            }
        }

        // 3. For each final block, compute x-center
        for (auto &blk : stack) {
            double sumX = 0.0;
            for (int i = blk.startIdx; i <= blk.endIdx; ++i) {
                sumX += x[i];
            }
            double len = static_cast<double>(blk.endIdx - blk.startIdx + 1);
            double center = sumX / len;
            // reuse 'sum' as x-center (or add a new field if you prefer)
            blk.sum = center;
        }

        // 4. Interpolate each final block's points
        std::vector<double> result(n, 0.0);

        for (size_t b = 0; b < stack.size(); ++b) {
            auto &curBlk  = stack[b];
            double curAvg = curBlk.avg;
            double curXc  = curBlk.sum; // x-center

            double nextAvg = curAvg;
            double nextXc  = curXc;
            if (b < stack.size() - 1) {
                auto &nb = stack[b+1];
                nextAvg = nb.avg;
                nextXc  = nb.sum;
            }

            for (int i = curBlk.startIdx; i <= curBlk.endIdx; ++i) {
                if (b == stack.size()-1 || std::fabs(nextXc - curXc) < 1e-15) {
                    // Last block or degenerate => constant
                    result[i] = curAvg;
                } else {
                    double t = (x[i] - curXc) / (nextXc - curXc);
                    if (t < 0.0) t = 0.0;
                    if (t > 1.0) t = 1.0;
                    result[i] = curAvg + t*(nextAvg - curAvg);
                }
                result[i] = max( min(result[i], max_y), min_y);
                assert((i < 1) || ( result[i-1] <= result[i] ));
            }
        }

        return result;
    }
};



/**
 * Inherits IsotonicRegression and
 * adds any specialized pre-processing or "PEP" logic.
 */
class IsotonicPEP : public IsotonicRegression
{
public:
    IsotonicPEP() = default;
    virtual ~IsotonicPEP() = default;

    std::vector<double> pepRegression(const std::vector<double>& y) const
    {
        std::vector<double> result = pavaNonDecreasingRanged(y, 0.0, 1.0);
        return result;
    }

    // Variant that produces a piecewise-linear function in logit space:
    std::vector<double> pepRegression(
        const std::vector<double>& y,
        const std::vector<double>& x
    ) const
    {
        if (x.size() != y.size()) {
            throw std::invalid_argument("pepRegressionWithInterpolation: x.size() != y.size()");
        }
        std::vector<double> result = pavaNonDecreasingInterpolationRanged(x, y, 0.0, 1.0);
        return result;
    }

    double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
        double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
        return interp_pep;
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
        std::vector<double> pep_iso = pepRegression(raw_pep);
        return pep_iso;
    }

    std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
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
        std::vector<double> pep_iso = pepRegression(raw_pep, scores);
        return pep_iso;
    }

    std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {}) {

        double epsilon=1e-20;
        auto is_dec = is_decoy;
        is_dec.insert(is_dec.begin(),0.5);

        std::vector<double> decoy_rate;
        if (!scores.empty()) {
            auto sc = scores;
            sc.insert(sc.begin(),sc[0]);
            decoy_rate = pepRegression(is_dec, sc);    
        } else {
            decoy_rate = pepRegression(is_dec);
        }
        decoy_rate.erase(decoy_rate.begin());

        std::vector<double> pep_iso;

        for (auto& dp : decoy_rate) {
            if (dp > 1. - epsilon)
                dp = 1. - epsilon;
            double pep = dp / ( 1 - dp );
            if (pep > 1.)
                pep = 1.;
            pep_iso.push_back(pep);
        }
        return pep_iso;
    }
   

protected:
    std::vector<double> qs;
    std::vector<double> pep_iso;

};
#endif /* ISOTONICPEP_H_ */
