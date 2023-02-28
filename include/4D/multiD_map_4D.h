#pragma once

#include <iterator>
#include <limits>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstring>
#include <math.h>
#include "rtree.h"
#include "rtree_objects.h"
#include "util.h"
#include <spatialindex/SpatialIndex.h>
#include "dkm/include/dkm.hpp"
#include <time.h>
#include "grid_file_builder_4D.h"

extern uint64_t arr_duration;
extern uint64_t arr_cnt;
std::vector<float> onerates;
namespace mdm
{
  template <class KeyType>
  class MultiDMap4D
  {
  public:
    MultiDMap4D(std::vector<std::vector<KeyType>> data, size_t num_radix_bits = 18, size_t cell_size = 100, size_t rs_group_size = 200, size_t exp = 0)
        : num_radix_bits_(num_radix_bits),
          rs_group_size_(rs_group_size),
          cell_size_(cell_size),
          exp_(exp)
    {
      size_t n = data.size();
      key_num_ = n;
      min_keys_ = util::GetMinKey<KeyType>(data, 4);
      max_keys_ = util::GetMaxKey<KeyType>(data, 4);
      group_size_ = ceil(pow(n / cell_size, 0.75) * cell_size);
      int cell_num_per_group_ = ceil(n / 1.0 / group_size_) + 1;
      data_.resize(n, std::vector<uint64_t>(4));
      ratio_ = pow(10, exp);
      for (int i = 0; i < n; i++)
      {
        data_[i][0] = round(data[i][0] * ratio_);
        data_[i][1] = round(data[i][1] * ratio_);
        data_[i][2] = round(data[i][2] * ratio_);
        data_[i][3] = round(data[i][3] * ratio_);
      }
      data_arr_.resize(cell_num_per_group_);
      for (int i = 0; i < cell_num_per_group_; i++)
      {
        data_arr_[i].resize(cell_num_per_group_);
        for (int j = 0; j < cell_num_per_group_; j++)
        {
          data_arr_[i][j].resize(cell_num_per_group_);
          for (int k = 0; k < cell_num_per_group_; k++)
          {
            data_arr_[i][j][k].resize(cell_num_per_group_);
          }
        }
      }
      std::sort(data_.begin(), data_.end(), [](const std::vector<uint64_t> &a, const std::vector<uint64_t> &b)
                { return a[0] < b[0]; });
      gf4D::Builder<uint64_t> gfb4d(data_.front()[0], data_.back()[0], group_size_, rs_group_size_, num_radix_bits_, cell_size_);
      for (size_t i = 0; i < n; i++)
      {
        if ((i + 1) % group_size_ == 0 || i == n - 1)
        {
          gfb4d.AddKeyToSpline(data_[i][0]);
        }
        vector<uint64_t> tmp = {data_[i][1], data_[i][2], data_[i][3]};
        gfb4d.AddKey(tmp);
      }
      gf4d_ = gfb4d.Finalize();
      float avgonerate = 0;
      for (int i = 0; i < onerates.size(); i++)
      {
        avgonerate += onerates[i];
      }
      avgonerate /= onerates.size();
      for (int i = 0; i < n; i++)
      {
        std::vector<uint_fast32_t> rtindex = gf4d_.GetCoord(data_[i]);
        data_arr_[rtindex[0]][rtindex[1]][rtindex[2]][rtindex[3]].push_back(data_[i]);
      }
      int cnt = 0;
      double err = 0;
      for (int i = 0; i < cell_num_per_group_; i++)
      {
        for (int j = 0; j < cell_num_per_group_; j++)
        {
          for (int k = 0; k < cell_num_per_group_; k++)
          {
            for (int l = 0; l < cell_num_per_group_; l++)
            {
              int size = data_arr_[i][j][k][l].size();
              cnt += size;
              err += abs(size - int(cell_size_));
            }
          }
        }
      }
      err /= pow(cell_num_per_group_, 4);
      std::cout << "size[MB]: " << (gf4d_.GetSizeInByte() / 1000) / 1000.0 << std::endl;
      std::cout << "max_prefix_all: " << max_prefix_all << std::endl;
      std::cout << "avgonerate: " << avgonerate << std::endl;
      std::cout << "accuracy: " << (cnt - gf::error_cnt) / 1.0 / cnt << std::endl;
      std::cout << "item_cnt " << cnt << std::endl;
      std::cout << "divide_err " << err << std::endl;
    }
    uint_fast64_t find(std::vector<KeyType> key)
    {
      std::vector<uint64_t> key_(4);
      std::vector<std::vector<uint64_t>> answ;
      key_[0] = round(key[0] * ratio_);
      key_[1] = round(key[1] * ratio_);
      key_[2] = round(key[2] * ratio_);
      key_[3] = round(key[3] * ratio_);
      std::vector<uint_fast32_t> rtindex = gf4d_.GetCoord(key_);
      auto arr_start = std::chrono::high_resolution_clock::now();
      for (std::vector<uint64_t> item : data_arr_[rtindex[0]][rtindex[1]][rtindex[2]][rtindex[3]])
      {
        if (item[0] == key_[0] && item[1] == key_[1] && item[2] == key_[2] && item[3] == key_[3])
        {
          answ.push_back(item);
        }
      }
      if (answ.size() == 0)
      {
        cerr << "can't find the key" << endl;
      }
      auto arr_end = std::chrono::high_resolution_clock::now();
      arr_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(arr_end - arr_start).count();
      arr_cnt += 1;
      total_cnt += answ.size();
      return answ[0][0];
    }

    void insert(std::vector<KeyType> key)
    {
      std::vector<uint64_t> key_(4);
      key_[0] = round(key[0] * ratio_);
      key_[1] = round(key[1] * ratio_);
      key_[2] = round(key[2] * ratio_);
      key_[3] = round(key[3] * ratio_);
      std::vector<uint_fast32_t> rtindex = gf4d_.GetCoord(key_);
      data_arr_[rtindex[0]][rtindex[1]][rtindex[2]][rtindex[3]].push_back(key_);
    }

    std::vector<std::vector<uint64_t>> range_search(std::vector<KeyType> begin, std::vector<KeyType> end)
    {
      std::vector<uint64_t> begin_(4);
      begin_[0] = round(begin[0] * ratio_);
      begin_[1] = round(begin[1] * ratio_);
      begin_[2] = round(begin[2] * ratio_);
      begin_[3] = round(begin[3] * ratio_);
      std::vector<uint64_t> end_(4);
      end_[0] = round(end[0] * ratio_);
      end_[1] = round(end[1] * ratio_);
      end_[2] = round(end[2] * ratio_);
      end_[3] = round(end[3] * ratio_);
      std::vector<std::vector<uint64_t>> total_answ;

      uint_fast32_t bx = gf4d_.GetCoordX(begin_[0]);
      uint_fast32_t ex = gf4d_.GetCoordX(end_[0]);
      for (int i = bx; i <= ex; i++)
      {
        uint_fast32_t by = gf4d_.GetCoordY(i, begin_[1]);
        uint_fast32_t ey = gf4d_.GetCoordY(i, end_[1]);
        for (int j = by; j <= ey; j++)
        {
          uint_fast32_t bz = gf4d_.GetCoordZ(i, j, begin_[2]);
          uint_fast32_t ez = gf4d_.GetCoordZ(i, j, end_[2]);
          for (int k = bz; k <= ez; k++)
          {
            uint_fast32_t bh = gf4d_.GetCoordH(i, j, k, begin_[3]);
            uint_fast32_t eh = gf4d_.GetCoordH(i, j, k, end_[3]);
            for (int l = bh; l <= eh; l++)
            {
              std::vector<std::vector<uint64_t>> answ;
              auto arr_start = std::chrono::high_resolution_clock::now();
              for (std::vector<uint64_t> item : data_arr_[i][j][k][l])
              {
                if (isPointInside(item, begin_, end_))
                {
                  answ.push_back(item);
                }
              }
              auto arr_end = std::chrono::high_resolution_clock::now();
              arr_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(arr_end - arr_start).count();
              arr_cnt += 1;
              total_answ.insert(total_answ.end(), answ.begin(), answ.end());
            }
          }
        }
      }
      // printf("%d\n",total_answ.size());
      total_cnt += total_answ.size();
      return total_answ;
    }

    std::vector<std::vector<uint64_t>> knn(std::vector<KeyType> key, int k_num)
    {
      float lenx = (max_keys_[0] - min_keys_[0]) * pow(k_num / 1.0 / key_num_, 0.25);
      float leny = (max_keys_[1] - min_keys_[1]) * pow(k_num / 1.0 / key_num_, 0.25);
      float lenz = (max_keys_[2] - min_keys_[2]) * pow(k_num / 1.0 / key_num_, 0.25);
      float lenh = (max_keys_[3] - min_keys_[3]) * pow(k_num / 1.0 / key_num_, 0.25);
      std::vector<KeyType> begin_(4);
      std::vector<KeyType> end_(4);
      std::vector<std::vector<uint64_t>> answ;
      while (true)
      {
        knn_query_cnt++;
        begin_[0] = max(float(0), key[0] - lenx / 2);
        begin_[1] = max(float(0), key[1] - leny / 2);
        begin_[2] = max(float(0), key[2] - lenz / 2);
        begin_[3] = max(float(0), key[3] - lenh / 2);
        end_[0] = key[0] + lenx / 2;
        end_[1] = key[1] + leny / 2;
        end_[2] = key[2] + lenz / 2;
        end_[3] = key[3] + lenh / 2;
        std::vector<std::vector<uint64_t>> answ = range_search(begin_, end_);
        if (answ.size() >= k_num)
        {
          break;
        }
        lenx *= pow(2, 0.25);
        leny *= pow(2, 0.25);
        lenz *= pow(2, 0.25);
        lenh *= pow(2, 0.25);
      }
      return answ;
    }

  private:
    bool isPointInside(std::vector<uint64_t> point, std::vector<uint64_t> begin, std::vector<uint64_t> end)
    {
      return point[0] >= begin[0] && point[0] <= end[0] && point[1] >= begin[1] && point[1] <= end[1] && point[2] >= begin[2] && point[2] <= end[2] && point[3] >= begin[3] && point[3] <= end[3];
    }
    std::vector<std::vector<uint64_t>> data_;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<uint64_t>>>>>> data_arr_;
    std::vector<float> min_keys_;
    std::vector<float> max_keys_;
    gf4D::GridFile<uint64_t> gf4d_;
    const size_t num_radix_bits_;
    size_t group_size_;
    size_t rs_group_size_;
    size_t cell_size_;
    size_t exp_;
    size_t ratio_;
    size_t key_num_;
  };
} // namespace mdm