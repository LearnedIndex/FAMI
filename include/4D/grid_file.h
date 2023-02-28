#pragma once
#include "libmorton/morton.h"
#include "radix_spline.h"
#include "common.h"
#include <vector>
#include <iostream>
extern uint64_t spline_duration;
namespace gf
{
  template <class KeyType>
  class GridFile
  {
  public:
    GridFile() = default;
// GridFile(std::vector<KeyType> min_keys, std::vector<KeyType> max_keys, int dim, double grid_size)
//     : min_keys_(min_keys), max_keys_(max_keys), dim_(dim), grid_size_(grid_size) {}
#ifndef USE_BITMAP
    GridFile(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t num_shift_bits, size_t group_size, size_t cell_size,
             std::vector<uint32_t> radix_table,
             std::vector<KeyType> spline_points,
             std::vector<rs::RadixSpline<KeyType>> rss)
        //  std::vector<MLP> mlps)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(num_shift_bits),
          group_size_(group_size),
          cell_size_(cell_size),
          dim_(2),
          radix_table_(std::move(radix_table)),
          spline_points_(std::move(spline_points)),
          rss_(std::move(rss))
    {
    }
#else
    GridFile(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t num_shift_bits, size_t group_size, size_t cell_size,
             std::bitset<BITMAP_SIZE> radix_bitmap,
             std::vector<uint32_t> radix_table,
             std::vector<KeyType> spline_points,
             std::vector<rs::RadixSpline<KeyType>> rss)
        //  std::vector<MLP> mlps)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(num_shift_bits),
          group_size_(group_size),
          cell_size_(cell_size),
          dim_(2),
          radix_bitmap_(std::move(radix_bitmap)),
          radix_table_(std::move(radix_table)),
          spline_points_(std::move(spline_points)),
          rss_(std::move(rss))
    {
    }
#endif

    // mlps_(std::move(mlps)) {}

    // Returns the size in bytes.
    size_t GetSizeInByte() const
    {
      size_t rs_size = 0;
      int cnt = 0;
      for (rs::RadixSpline<KeyType> rs : rss_)
      {
        rs_size += rs.GetSize();
        cnt++;
      }
      return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) + spline_points_.size() * sizeof(KeyType) + rs_size;
    }

    std::vector<uint_fast32_t> GetCoord(std::vector<KeyType> key)
    {
      if (key[0] < min_key_)
        return {0, 0};
      if (key[0] > max_key_)
        return {spline_points_.size() - 1, 0};

      // auto spline_start = std::chrono::high_resolution_clock::now();
      const size_t index = GetSplineSegment(key[0]);
      // auto spline_end = std::chrono::high_resolution_clock::now();
      // spline_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(spline_end - spline_start).count();

      // torch::Tensor input = torch::tensor({float(key[1])});
      // input = torch::unsqueeze(input, 1);
      // // std::cout << input << std::endl;
      // double ret = mlps_[index].forward(input).item().toDouble();
      // if (ret < 0)
      //   return {index, 0};
      // if (ret >= (group_size_ - 1) / cell_size_)
      //   return {index, (group_size_ - 1) / cell_size_};
      // return {index, ret};

      // double ret = rss_[index].GetEstimatedPosition(key[1]);
      // if (ret < 0)
      //   return {index, 0};
      // if (ret >= group_size_)
      //   return {index, (group_size_ - 1) / cell_size_};
      // return {index, ret / cell_size_};

      double ret = rss_[index].GetEstimatedPosition(key[1]);
      return {index, ret};
    }

    uint_fast32_t GetCoordX(KeyType keyx)
    {
      if (keyx < min_key_)
        return 0;
      if (keyx > max_key_)
        return spline_points_.size() - 1;
      return GetSplineSegment(keyx);
    }

    uint_fast32_t GetCoordY(uint_fast32_t x, KeyType keyy)
    {
      double ret = rss_[x].GetEstimatedPosition(keyy);
      return ret;
    }

    float GetRadixTableAvgDis()
    {
      float res = 0;
      for (int i = 1; i < radix_table_.size(); i++)
      {
        res += radix_table_[i] - radix_table_[i - 1];
      }
      return res / radix_table_.size();
    }

  private:
    // std::vector<uint_fast32_t> GetCoord(std::vector<KeyType> key)
    // {
    //     std::vector<uint_fast32_t> coord(dim_);
    //     for (int i = 0; i < dim_; i++)
    //     {
    //         coord[i] = (key[i] - min_keys_[i]) / grid_size_;
    //     }
    //     return coord;
    // }

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
    // std::vector<KeyType> min_keys_;
    // std::vector<KeyType> max_keys_;
    // double grid_size_;

    int dim_;

    KeyType min_key_;
    KeyType max_key_;
    size_t num_radix_bits_;
    size_t num_shift_bits_;
    size_t group_size_;
    size_t cell_size_;
#ifdef USE_BITMAP
    std::bitset<BITMAP_SIZE> radix_bitmap_;
#endif
    std::vector<uint32_t> radix_table_;
    std::vector<KeyType> spline_points_;
    std::vector<rs::RadixSpline<KeyType>> rss_;

  };
}