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
     * on the input 'values'.  The 'counts' array indicates how many times each 'value'
     * is repeated. If omitted, each point is assumed to have count=1.
     *
     * Returns an expanded vector with length = sum(counts).
     * Each final "block" is repeated block.count times in the output.
     *
     * Typically used in 1D for y-values only.
     */
    virtual std::vector<double> pavaNonDecreasing(
        const std::vector<double>& values,
        const std::vector<int>& counts
    ) const
    {
        if (values.size() != counts.size()) {
            throw std::invalid_argument("pavaNonDecreasing: values.size() != counts.size()");
        }
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
            double sum_i = values[i] * counts[i];
            MergeBlock newBlock { sum_i, counts[i], values[i] };
            stack.push_back(newBlock);

            // 2. Merge while there's a violation of non-decreasing
            while (stack.size() > 1) {
                auto &top    = stack.back();
                auto &secTop = stack[stack.size() - 2];

                if (secTop.avg > top.avg) {
                    double mergedSum   = secTop.sum + top.sum;
                    int    mergedCount = secTop.count + top.count;
                    double mergedAvg   = mergedSum / mergedCount;
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
            for (int c = 0; c < b.count; ++c) {
                result.push_back(b.avg);
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
    virtual std::vector<double> pavaNonDecreasingInterpolation(
        const std::vector<double>& x,
        const std::vector<double>& y
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
        std::vector<ExtBlock> stack;
        stack.reserve(n);

        for (auto &b : blocks) {
            stack.push_back(b);

            // Merge while there's a violation
            while (stack.size() > 1) {
                auto &top    = stack.back();
                auto &secTop = stack[stack.size()-2];
                if (secTop.avg > top.avg) {
                    double newSum   = secTop.sum + top.sum;
                    int    newCount = secTop.count + top.count;
                    double newAvg   = newSum / newCount;
                    int    sIdx     = secTop.startIdx;
                    int    eIdx     = top.endIdx;
                    stack.pop_back();
                    stack.pop_back();
                    stack.push_back({
                        newSum, newCount, newAvg, sIdx, eIdx
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
            }
        }

        return result;
    }
};


/**
 * Derived class that applies logistic transformations
 * around the base IsotonicRegression’s PAVA routines.
 *
 * - We have clampProbability(), logistic(), logit()
 * - logisticIsotonicRegression(...) => stepwise-constant in logit space
 * - logisticIsotonicInterpolation(...) => piecewise-linear in logit space
 */
class LogisticIsotonicRegression : public IsotonicRegression
{
protected:
    static double clampProbability(double p, double eps = 1e-12)
    {
        if (p < eps)       return eps;
        if (p > 1.0 - eps) return 1.0 - eps;
        return p;
    }

    static double logistic(double x)
    {
        // 1 / (1 + e^-x)
        return 1.0 / (1.0 + std::exp(-x));
    }

    static double logit(double p)
    {
        p = clampProbability(p);
        return std::log(p / (1.0 - p));
    }

public:
    LogisticIsotonicRegression() = default;
    virtual ~LogisticIsotonicRegression() = default;

    /**
     * logisticIsotonicRegression:
     *   - Input y-values (which might be outside [0,1])
     *   - Optionally "clip" or "merge" them externally if needed
     *   - Convert each y_i to logit, run pavaNonDecreasing in logit-space,
     *     then convert back with logistic().
     *
     * If 'counts' is empty, we assume each y has count=1.
     */
    virtual std::vector<double> logisticIsotonicRegression(
        const std::vector<double>& y,
        const std::vector<int>& counts = {}
    ) const
    {
        const int n = static_cast<int>(y.size());
        if (n == 0) {
            return {};
        }

        // If no counts given, assume 1 for each
        std::vector<int> c;
        if (counts.empty()) {
            c.resize(n, 1);
        } else {
            c = counts;
        }

        // 1. Convert y => logit (with clamp)
        std::vector<double> logitVals(n);
        for (int i = 0; i < n; ++i) {
            double p = clampProbability(y[i]);
            logitVals[i] = logit(p);
        }

        // 2. PAVA in logit space
        std::vector<double> mergedLogits = pavaNonDecreasing(logitVals, c);
        // mergedLogits has size = sum(c).

        // 3. Convert back => logistic
        std::vector<double> result(mergedLogits.size());
        for (size_t i = 0; i < mergedLogits.size(); ++i) {
            result[i] = logistic(mergedLogits[i]);
        }

        return result;
    }


    /**
     * logisticIsotonicInterpolation:
     *   - Similar idea, but we do the "interpolation" variant in logit space.
     *   - We rely on pavaNonDecreasingInterpolation(...) in the *base class*,
     *     but we pass it the logit-transformed values. Then we re-transform back
     *     with logistic().
     *
     *   Here we need x as well, for the interpolation.
     *
     *   In actual usage, you might have more steps, e.g. merging blocks until
     *   they're in [0,1]. But let's keep it simpler here.
     */
    virtual std::vector<double> logisticIsotonicInterpolation(
        const std::vector<double>& x,
        const std::vector<double>& y
    ) const
    {
        if (x.size() != y.size()) {
            throw std::invalid_argument("logisticIsotonicInterpolation: x.size() != y.size()");
        }
        const int n = static_cast<int>(y.size());
        if (n == 0) {
            return {};
        }

        // 1. logit-transform y
        std::vector<double> logitVals(n);
        for (int i = 0; i < n; ++i) {
            double p = clampProbability(y[i]);
            logitVals[i] = logit(p);
        }

        // 2. Use the base interpolation method, which merges in real space
        //    but we treat "logitVals" as the data to be made isotonic.
        //    So the "non-decreasing" constraint is in logit space.
        //    We pass x, logitVals to pavaNonDecreasingInterpolation.
        //
        //    This will produce a piecewise-linear function in logit space
        //    across the x-centers.
        std::vector<double> logitInterp = pavaNonDecreasingInterpolation(x, logitVals);

        // 3. Convert back => logistic
        std::vector<double> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = logistic(logitInterp[i]);
        }
        return result;
    }
};


/**
 * Inherits from LogisticIsotonicRegression and
 * adds any specialized pre-processing or "PEP" logic.
 *
 * For example:
 *   - Merging consecutive y until each block average is in [0,1]
 *   - Then calling logisticIsotonicRegression(...) from the parent
 *     or logisticIsotonicInterpolation(...), etc.
 */
class IsotonicPEP : public LogisticIsotonicRegression
{
public:
    IsotonicPEP() = default;
    virtual ~IsotonicPEP() = default;

    /**
     * Example: "pepRegression" that merges consecutive points
     * so that each block has an average in [0,1], then does
     * logistic isotonic regression.
     */
    std::vector<double> pepRegression(const std::vector<double>& y) const
    {
        // 1. Merge consecutive points so each block's average is in [0,1].
        std::vector<double> mergedY = createBlocksInUnitIntervalAndUnfold(y);

        // 2. Now call parent's logistic method on the merged data.
        //    The data is length(mergedY) now. Each point has count=1.
        std::vector<double> result = logisticIsotonicRegression(mergedY);
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
        // 1) Merge consecutive y so each block average is in [0,1]
        std::vector<double> mergedY = createBlocksInUnitIntervalAndUnfold(y);

        // 2) Interpolate in logit space. This uses x for the piecewise-linear stepping.
        //    'logisticIsotonicInterpolation' merges in logit space, then
        //    does linear interpolation between block centers (in x).
        std::vector<double> result = logisticIsotonicInterpolation(x, mergedY);
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

protected:
    /**
     * Merge consecutive y[i] so that each partial block average is in [0,1].
     * Then "unfold" back into a single vector. (This is a simplified approach.)
     *
     * For example, if you have data:
     *   y = [-0.2, 0.1, 0.05, 1.2, 0.5]
     * you might accumulate: 
     *   block1 => [-0.2, 0.1, 0.05] => sum= -0.05 => avg= -0.0167 => still out of [0,1]
     *   (just a conceptual example, you might keep going until it falls in range or clip).
     *
     * The approach shown: "accumulate until average is in [0,1], then finalize block".
     * If leftover is out, we clip it. Finally, we flatten blocks back into a single vector.
     */
    std::vector<double> createBlocksInUnitIntervalAndUnfold(const std::vector<double>& y) const
    {
        std::vector<std::vector<double>> blocks;
        blocks.reserve(y.size());

        double currentSum   = 0.0;
        int    currentCount = 0;

        std::vector<double> currentBlock;
        currentBlock.reserve(y.size()); 

        for (double val : y) {
            currentSum   += val;
            currentCount += 1;
            double avg = currentSum / currentCount;
            if (avg > 0.0 && avg < 1.0) {
                std::vector<double> currentBlock(currentCount, avg);
                // finalize
                blocks.push_back(currentBlock);
                // reset
                currentSum   = 0.0;
                currentCount = 0;
                currentBlock.clear();
            }
        }

        // leftover
        if (currentCount>0) {
            double avg = currentSum / currentCount;
            // If out of [0,1], clip
            if (avg < 0.0) avg = 1e-10;
            if (avg > 1.0) avg = 1.0 - 1e-10;

            // Instead of returning single avg, we might simply fill them with 'avg'
            // or just store them as is. You can choose your approach.
            // For simplicity, let's store them all as 'avg'.
            std::vector<double> clippedBlock(currentCount, avg);
            blocks.push_back(clippedBlock);
        }

        // "Unfold" all blocks into a single vector
        std::vector<double> result;
        for (auto &b : blocks) {
            result.insert(result.end(), b.begin(), b.end());
        }

        return result;
    }

protected:
    std::vector<double> qs;
    std::vector<double> pep_iso;

};
#endif /* ISOTONICPEP_H_ */
