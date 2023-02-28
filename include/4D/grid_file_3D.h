#pragma once
#include "libmorton/morton.h"
#include "radix_spline.h"
#include "common.h"
#include "grid_file.h"
#include <vector>
#include <iostream>
extern uint64_t spline_duration;
namespace gf3D
{
  template <class KeyType>
  class GridFile
  {
  public:
    GridFile() = default;
// GridFile(std::vector<KeyType> min_keys, std::vector<KeyType> max_keys, int dim, double grid_size)
//     : min_keys_(min_keys), max_keys_(max_keys), dim_(dim), grid_size_(grid_size) {}
#ifndef USE_BITMAP
    GridFile(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t num_shift_bits,
             std::vector<uint32_t> radix_table,
             std::vector<KeyType> spline_points,
             std::vector<gf::GridFile<uint64_t>> gfs)
        //  std::vector<MLP> mlps)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(num_shift_bits),
          radix_table_(std::move(radix_table)),
          spline_points_(std::move(spline_points)),
          gfs_(std::move(gfs))
    {
    }
#else
    GridFile(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t num_shift_bits,
             std::bitset<BITMAP_SIZE> radix_bitmap,
             std::vector<uint32_t> radix_table,
             std::vector<KeyType> spline_points,
             std::vector<gf::GridFile<uint64_t>> gfs)
        //  std::vector<MLP> mlps)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(num_shift_bits),
          radix_bitmap_(std::move(radix_bitmap)),
          radix_table_(std::move(radix_table)),
          spline_points_(std::move(spline_points)),
          gfs_(std::move(gfs))
    {
    }
#endif

    // Returns the size in bytes.
    size_t GetSizeInByte() const
    {
      size_t gf_size = 0;
      int cnt = 0;
      for (gf::GridFile<uint64_t> gf : gfs_)
      {
        gf_size += gf.GetSizeInByte();
        cnt++;
      }
      return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) + spline_points_.size() * sizeof(KeyType) + gf_size;
    }

    std::vector<uint_fast32_t> GetCoord(std::vector<KeyType> key)
    {
      if (key[0] < min_key_)
        return {0, 0, 0};
      if (key[0] > max_key_)
        return {spline_points_.size() - 1, 0, 0};

      // auto spline_start = std::chrono::high_resolution_clock::now();
      const size_t index = GetSplineSegment(key[0]);
      const vector<size_t> ret = gfs_[index].GetCoord({key[1],key[2]});
      return {index, ret[0],ret[1]};
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
      const size_t ret = gfs_[x].GetCoordX(keyy);
      return ret;
    }

    uint_fast32_t GetCoordZ(uint_fast32_t x, uint_fast32_t y, KeyType keyz)
    {
      const size_t ret = gfs_[x].GetCoordY(y,keyz);
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
    size_t num_radix_bits_;
    size_t num_shift_bits_;
#ifdef USE_BITMAP
    std::bitset<BITMAP_SIZE> radix_bitmap_;
#endif
    std::vector<uint32_t> radix_table_;
    std::vector<KeyType> spline_points_;
    std::vector<gf::GridFile<uint64_t>> gfs_;
  };
}