#pragma once

#include <cassert>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
// #include "MLP.h"
#include "common.h"
#include "grid_file_builder.h"
#include "grid_file_3D.h"
extern std::vector<float> onerates;
extern uint32_t max_prefix_all;
namespace gf3D
{
  int error_cnt = 0;
  template <class KeyType>
  class Builder
  {
  public:
    Builder(KeyType min_key, KeyType max_key, size_t group_size, size_t rs_group_size, size_t num_radix_bits = 18, size_t cell_size = 100)
        : min_key_(min_key),
          max_key_(max_key),
          group_size_(group_size),
          rs_group_size_(rs_group_size),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(GetNumShiftBits(max_key - min_key, num_radix_bits_)),
          cell_size_(cell_size),
          prev_prefix_(-1)
    {
      // Initialize radix table, needs to contain all prefixes up to the largest key + 1.
      const uint32_t max_prefix = (max_key - min_key) >> num_shift_bits_;
      max_prefix_ = max_prefix;
      max_prefix_all = std::max(max_prefix_all,max_prefix_);
#ifndef USE_BITMAP
      radix_table_.resize(max_prefix + 2, 0);
#endif
    }

    void AddKey(vector<KeyType> key)
    {
      group_points_.push_back(key);
      if (group_points_.size() == group_size_)
      {
        DoFitting();
        group_points_.clear();
      }
    }

    void AddKeyToSpline(KeyType key)
    {
      spline_points_.push_back(key);
      PossiblyAddKeyToRadixTable(key);
    }

    // Finalizes the construction and returns a read-only `GridFile`.
    gf3D::GridFile<KeyType> Finalize()
    {
      if (group_points_.size() != 0)
      {
        DoFitting();
        group_points_.clear();
      }
      // Maybe even size the radix based on max key right from the start
      FinalizeRadixTable();
#ifndef USE_BITMAP
      return GridFile<KeyType>(min_key_,
                               max_key_,
                               num_radix_bits_,
                               num_shift_bits_,
                               std::move(radix_table_),
                               std::move(spline_points_),
                               std::move(gfs_));
#else
      // printf("rs onerate: %lf\n", radix_bitmap_.count() / 1.0 / max_prefix_);
      onerates.push_back(radix_bitmap_.count() / 1.0 / max_prefix_);
      return gf3D::GridFile<KeyType>(min_key_,
                                     max_key_,
                                     num_radix_bits_,
                                     num_shift_bits_,
                                     std::move(radix_bitmap_),
                                     std::move(radix_table_),
                                     std::move(spline_points_),
                                     std::move(gfs_));
#endif
      //  std::move(mlps_));
    }

  private:
    // Returns the number of shift bits based on the `diff` between the largest and the smallest key.
    // KeyType == uint32_t.
    static size_t GetNumShiftBits(uint32_t diff, size_t num_radix_bits)
    {
      const uint32_t clz = __builtin_clz(diff);
      if ((32 - clz) < num_radix_bits)
        return 0;
      return 32 - num_radix_bits - clz;
    }
    // KeyType == uint64_t.
    static size_t GetNumShiftBits(uint64_t diff, size_t num_radix_bits)
    {
      const uint32_t clzl = __builtin_clzl(diff);
      if ((64 - clzl) < num_radix_bits)
        return 0;
      return 64 - num_radix_bits - clzl;
    }

    void PossiblyAddKeyToRadixTable(KeyType key)
    {
#ifndef USE_BITMAP
      const KeyType curr_prefix = (key - min_key_) >> num_shift_bits_;
      if (curr_prefix != prev_prefix_)
      {
        const uint32_t curr_index = spline_points_.size() - 1;
        for (KeyType prefix = prev_prefix_ + 1; prefix <= curr_prefix; ++prefix)
          radix_table_[prefix] = curr_index;
        prev_prefix_ = curr_prefix;
      }
#else
      const KeyType curr_prefix = (key - min_key_) >> num_shift_bits_;
      if (curr_prefix != prev_prefix_)
      {
        const uint32_t curr_index = spline_points_.size() - 1;
        radix_bitmap_.set(curr_prefix);
        radix_table_.push_back(curr_index);
        prev_prefix_ = curr_prefix;
      }
#endif
    }

    void FinalizeRadixTable()
    {
#ifndef USE_BITMAP
      ++prev_prefix_;
      const uint32_t num_spline_points = spline_points_.size();
      for (; prev_prefix_ < radix_table_.size(); ++prev_prefix_)
        radix_table_[prev_prefix_] = num_spline_points;
#else
      const uint32_t num_spline_points = spline_points_.size();
      // radix_bitmap_.set(radix_table_.size() - 1);
      radix_table_.push_back(num_spline_points - 1);
#endif
    }

    void DoFitting()
    {
      std::sort(group_points_.begin(), group_points_.end(), [](const std::vector<uint64_t> &a, const std::vector<uint64_t> &b)
                { return a[0] < b[0]; });
      int n = group_points_.size();
      int group_size_2d_ = ceil(group_size_ / pow(group_size_ / cell_size_, 1.0 / 2));
      gf::Builder<uint64_t> gfb(group_points_.front()[0], group_points_.back()[0], group_size_2d_, rs_group_size_, num_radix_bits_, cell_size_);
      for (size_t i = 0; i < n; i++)
      {
        if ((i + 1) % group_size_2d_ == 0 || i == n - 1)
        {
          gfb.AddKeyToSpline(group_points_[i][0]);
        }
        gfb.AddKey(group_points_[i][1]);
      }
      gf::GridFile<uint64_t> gf_ = gfb.Finalize();
      gfs_.push_back(gf_);
    }

    const KeyType min_key_;
    const KeyType max_key_;
    const size_t num_radix_bits_;
    const size_t num_shift_bits_;
    std::vector<uint32_t> radix_table_;
#ifdef USE_BITMAP
    std::bitset<BITMAP_SIZE> radix_bitmap_;
#endif
    std::vector<KeyType> spline_points_;
    std::vector<gf::GridFile<uint64_t>> gfs_;

    // std::vector<MLP> mlps_;

    KeyType prev_prefix_;

    uint32_t max_prefix_;

    size_t rs_group_size_;

    size_t cell_size_;

    // Number of CDF points in one group
    size_t group_size_;
    // vector of a CDF points group
    std::vector<std::vector<KeyType>> group_points_;
  };
} // namespace gf