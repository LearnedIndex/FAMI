#pragma once

#include <algorithm>
#include <vector>
#include <iostream>
#include <bitset>
#include "common.h"
// #define BITMAP_SIZE 512
extern uint64_t spline_duration;
namespace rs
{

    // Approximates a cumulative distribution function (CDF) using spline interpolation.
    template <class KeyType>
    class RadixSpline
    {
    public:
        RadixSpline() = default;
#ifndef USE_BITMAP
        RadixSpline(KeyType min_key,
                    KeyType max_key,
                    size_t num_keys,
                    size_t num_radix_bits,
                    size_t num_shift_bits,
                    size_t group_size,
                    std::vector<uint32_t> radix_table,
                    std::vector<KeyType> spline_points,
                    std::vector<FitParam> fit_params)
            : min_key_(min_key),
              max_key_(max_key),
              num_keys_(num_keys),
              num_radix_bits_(num_radix_bits),
              num_shift_bits_(num_shift_bits),
              group_size_(group_size),
              radix_table_(std::move(radix_table)),
              spline_points_(std::move(spline_points)),
              fit_params_(std::move(fit_params))
        {
        }
#else
        RadixSpline(KeyType min_key,
                    KeyType max_key,
                    size_t num_keys,
                    size_t num_radix_bits,
                    size_t num_shift_bits,
                    size_t group_size,
                    std::bitset<BITMAP_SIZE> radix_bitmap,
                    std::vector<uint32_t> radix_table,
                    std::vector<KeyType> spline_points,
                    std::vector<FitParam> fit_params)
            : min_key_(min_key),
              max_key_(max_key),
              num_keys_(num_keys),
              num_radix_bits_(num_radix_bits),
              num_shift_bits_(num_shift_bits),
              group_size_(group_size),
              radix_bitmap_(std::move(radix_bitmap)),
              radix_table_(std::move(radix_table)),
              spline_points_(std::move(spline_points)),
              fit_params_(std::move(fit_params))
        {
        }
#endif
        // Returns the estimated position of `key`.
        double GetEstimatedPosition(const KeyType key) const
        {
            // Truncate to data boundaries.
            if (key <= min_key_)
                return 0;
            if (key >= max_key_)
                return num_keys_ - 1;

            // Find spline segment with `key` ∈ (spline[index - 1], spline[index]].
            // auto spline_start = std::chrono::high_resolution_clock::now();
            const size_t index = GetSplineSegment(key);
            // auto spline_end = std::chrono::high_resolution_clock::now();
            // spline_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(spline_end - spline_start).count();

            const FitParam param = fit_params_[index];
            /* Compute estimate value */
            const double slope = param.slope;
            const double intercept = param.intercept;
            // Interpolate.
            double ret = key * slope + intercept;
            // double ret = std::fma(key_diff, slope, y);
            ret = std::max(group_size_ * index - 1.0, ret);
            if (ret < 0)
                return 0;
            if (ret >= num_keys_)
                return num_keys_ - 1;
            return ret;
        }

        // Returns the size in bytes.
        size_t GetSize() const
        {
            return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) + spline_points_.size() * sizeof(KeyType) + fit_params_.size() * sizeof(FitParam);
        }

        size_t GetSplineNum() const
        {
            return spline_points_.size();
        }

    private:
        // Returns the index of the spline point that marks the end of the spline segment that contains the `key`:
        // `key` ∈ (spline[index - 1], spline[index]]
        size_t GetSplineSegment(const KeyType key) const
        {
            // Narrow search range using radix table.
#ifndef USE_BITMAP
            const KeyType prefix = (key - min_key_) >> num_shift_bits_;
            assert(prefix + 1 < radix_table_.size());
            const uint32_t begin = radix_table_[prefix];
            const uint32_t end = radix_table_[prefix + 1];
#else
            KeyType prefix = (key - min_key_) >> num_shift_bits_;
            const uint32_t begin = radix_table_[(radix_bitmap_ << BITMAP_SIZE - prefix).count()];
            const uint32_t end = radix_table_[(radix_bitmap_ << BITMAP_SIZE - prefix - 1).count()];
#endif
            // knn_query_cnt += end - begin;
            if (end - begin < 32)
            {
                // Do linear search over narrowed range.
                uint32_t current = begin;
                while (spline_points_[current] < key)
                    ++current;
                return current;
            }

            // Do binary search over narrowed range.
            const auto lb = std::lower_bound(spline_points_.begin() + begin,
                                             spline_points_.begin() + end,
                                             key);
            return std::distance(spline_points_.begin(), lb);
        }

        KeyType min_key_;
        KeyType max_key_;
        size_t num_keys_;
        size_t num_radix_bits_;
        size_t num_shift_bits_;
        size_t group_size_;

#ifdef USE_BITMAP
        std::bitset<BITMAP_SIZE> radix_bitmap_;
#endif
        std::vector<uint32_t> radix_table_;
        std::vector<KeyType> spline_points_;
        std::vector<FitParam> fit_params_;
    };

} // namespace rs
