#pragma once

#include <iterator>
#include <limits>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstring>
#include <math.h>
#include "grid_file_builder.h"
#include "grid_file.h"
#include "rtree.h"
#include "rtree_objects.h"
#include "util.h"
#include <spatialindex/SpatialIndex.h>
#include "dkm/include/dkm.hpp"
#include <time.h>

extern uint64_t rtree_duration;
extern uint64_t rtree_cnt;
std::vector<float> onerates;
namespace mdm
{
  class MyVisitor : public SpatialIndex::IVisitor
  {
  public:
    void visitNode(const SpatialIndex::INode & /* n */) override {}

    void visitData(const SpatialIndex::IData &d) override
    {
      // std::cout << d.getIdentifier() << std::endl;
      rstar_res.push_back(d.getIdentifier());
      // the ID of this data entry is an answer to the query. I will just print it to stdout.
    }

    void visitData(std::vector<const SpatialIndex::IData *> & /* v */) override {}
  };

  template <class KeyType>
  class MultiDMap
  {
  public:
    MultiDMap(std::vector<std::vector<KeyType>> data, size_t num_radix_bits = 18, size_t cell_size = 100, size_t rs_group_size = 200, size_t exp = 0)
        : num_radix_bits_(num_radix_bits),
          group_size_(std::sqrt(data.size() * cell_size)),
          rs_group_size_(rs_group_size),
          cell_size_(cell_size),
          exp_(exp)
    {
      size_t n = data.size();
      key_num_ = n;
      min_keys_ = util::GetMinKey<KeyType>(data, 2);
      max_keys_ = util::GetMaxKey<KeyType>(data, 2);
      group_size_ = sqrt(n * cell_size_);
      data_.resize(n, std::vector<uint64_t>(2));
      rtree_size_ = group_size_ / cell_size_ + 1;
      rtrees_.resize(rtree_size_, std::vector<RTree<uint64_t>>(rtree_size_));
      for (int i = 0; i < rtree_size_; i++)
      {
        for (int j = 0; j < rtree_size_; j++)
        {
          rtrees_[i][j] = RTree<uint64_t>(40); // 64M需改为30
        }
      }
      ratio_ = pow(10, exp);
      for (int i = 0; i < n; i++)
      {
        data_[i][0] = round(data[i][0] * ratio_);
        data_[i][1] = round(data[i][1] * ratio_);
      }
      std::sort(data_.begin(), data_.end(), [](const std::vector<uint64_t> &a, const std::vector<uint64_t> &b)
                { return a[0] < b[0]; });
      gf::Builder<uint64_t> gfb(data_.front()[0], data_.back()[0], group_size_, rs_group_size_, num_radix_bits_, cell_size_);
      for (size_t i = 0; i < n; i++)
      {
        if ((i + 1) % group_size_ == 0 || i == n - 1)
        {
          gfb.AddKeyToSpline(data_[i][0]);
        }
        gfb.AddKey(data_[i][1]);
        pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(data_[i][0], data_[i][1])));
      }
      gf_ = gfb.Finalize();
      float avgonerate = 0;
      for (int i = 0; i < onerates.size(); i++)
      {
        avgonerate += onerates[i];
      }
      avgonerate /= onerates.size();
      float radixavgdis = gf_.GetRadixTableAvgDis();
      for (int i = 0; i < pols_.size(); i++)
      {
        // std::cout << i << std::endl;
        std::vector<uint_fast32_t> rtindex = gf_.GetCoord({pols_[i].get_Pmin().get_X(), pols_[i].get_Pmin().get_Y()});

        rtrees_[rtindex[0]][rtindex[1]].insert_polygon(&pols_[i], &pols_[i]);

        // if (!rtrees_.count(index))
        // {
        //   // std::string baseName = "Rtrees/Rtree" + std::to_string(index);
        //   // SpatialIndex::IStorageManager *diskfile = SpatialIndex::StorageManager::createNewDiskStorageManager(baseName, 4096);
        //   // SpatialIndex::StorageManager::IBuffer *file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);

        //   SpatialIndex::IStorageManager *memoryfile = SpatialIndex::StorageManager::createNewMemoryStorageManager();
        //   SpatialIndex::StorageManager::IBuffer *file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*memoryfile, 10, false);
        //   indexIdentifier_.insert({index,0});
        //   SpatialIndex::ISpatialIndex *tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier_[index]);
        //   rtrees_.insert({index, tree});
        // }
        // std::string value = std::to_string(pols_[i].get_Pmin().get_X());
        // double p[2];
        // p[0] = pols_[i].get_Pmin().get_X();
        // p[1] = pols_[i].get_Pmin().get_Y();
        // SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
        // rtrees_.find(index)->second->insertData((uint32_t)(value.size() + 1), reinterpret_cast<const uint8_t *>(value.c_str()), r, pols_[i].get_Pmin().get_X());
      }
      int cnt = 0;
      double err = 0;
      for (int i = 0; i < rtree_size_; i++)
      {
        for (int j = 0; j < rtree_size_; j++)
        {
          cnt += rtrees_[i][j].get_items_cnt();
          err += abs(rtrees_[i][j].get_items_cnt() - int(cell_size_));
        }
      }
      // std::cout << rtrees_[0][0].show_values_JSON() << std::endl;
      err /= rtrees_.size() * rtrees_[0].size();
      std::cout << "max_prefix_all: " << max_prefix_all << std::endl;
      std::cout << "avgonerate: " << avgonerate << std::endl;
      std::cout << "size[MB]: " << (gf_.GetSizeInByte() / 1000.0) / 1000.0 << std::endl;
      std::cout << "accuracy: " << (cnt - gf::error_cnt) / 1.0 / cnt << std::endl;
      std::cout << "item_cnt " << cnt << std::endl;
      std::cout << "divide_err " << err << std::endl;
    }
    uint_fast64_t find(std::vector<KeyType> key)
    {
      std::vector<uint64_t> key_(2);
      key_[0] = round(key[0] * ratio_);
      key_[1] = round(key[1] * ratio_);
      // return gf_.GetGridIndex(key_);
      std::vector<uint_fast32_t> rtindex = gf_.GetCoord(key_);
      // return rtindex[0];
      std::vector<Polygon<uint64_t> *> answ;
      // auto rtree_start = std::chrono::high_resolution_clock::now();
      rtrees_[rtindex[0]][rtindex[1]].range_search(Polygon<uint64_t>(Point<uint64_t>(key_[0], key_[1])), answ);
      // auto rtree_end = std::chrono::high_resolution_clock::now();
      // rtree_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(rtree_end - rtree_start).count();
      // rtree_cnt += 1;
      total_cnt += answ.size();
      // if (answ[0]->get_Pmin().get_X() != key_[0] || answ[0]->get_Pmin().get_Y() != key_[1])
      // {
      //   printf("error result\n");
      //   exit(-1);
      // }
      return answ[0]->get_Pmin().get_X();

      // double p[2];
      // p[0] = key_[0];
      // p[1] = key_[1];
      // MyVisitor vis;
      // SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
      // rtrees_.find(rtindex)->second->intersectsWithQuery(r, vis);
      // // if(res!=key_[0])
      // // {
      // //   printf("error\n");
      // //   // exit(-1);
      // // }
      // return res;
    }

    void insert(std::vector<KeyType> key)
    {
      std::vector<uint64_t> key_(2);
      key_[0] = round(key[0] * ratio_);
      key_[1] = round(key[1] * ratio_);
      std::vector<uint_fast32_t> rtindex = gf_.GetCoord(key_);

      pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(key_[0], key_[1])));
      rtrees_[rtindex[0]][rtindex[1]].insert_polygon(&pols_.back(), &pols_.back());
    }

    std::vector<Polygon<uint64_t> *> range_search(std::vector<KeyType> begin, std::vector<KeyType> end)
    {
      std::vector<uint64_t> begin_(2);
      begin_[0] = round(begin[0] * ratio_);
      begin_[1] = round(begin[1] * ratio_);
      std::vector<uint64_t> end_(2);
      end_[0] = round(end[0] * ratio_);
      end_[1] = round(end[1] * ratio_);
      std::vector<Polygon<uint64_t> *> total_answ;

      // std::vector<uint_fast32_t> btindex = gf_.GetCoord(begin_);
      // std::vector<uint_fast32_t> etindex = gf_.GetCoord(end_);
      // for (int i = btindex[0]; i <= etindex[0]; i++)
      // {
      //   for (int j = btindex[1]; j <= etindex[1]; j++)
      //   {
      //     std::vector<Polygon<uint64_t> *> answ;
      //     rtrees_[i][j].range_search(Polygon<uint64_t>(Point<uint64_t>(begin_[0], begin_[1]), Point<uint64_t>(end_[0], end_[1])), answ);
      //     total_answ.insert(total_answ.end(), answ.begin(), answ.end());
      //   }
      // }

      uint_fast32_t bx = gf_.GetCoordX(begin_[0]);
      uint_fast32_t ex = gf_.GetCoordX(end_[0]);
      for (int i = bx; i <= ex; i++)
      {
        uint_fast32_t by = gf_.GetCoordY(i, begin_[1]);
        uint_fast32_t ey = gf_.GetCoordY(i, end_[1]);
#ifndef ONLY_GRID
        for (int j = by; j <= ey; j++)
        {
          std::vector<Polygon<uint64_t> *> answ;
          // auto rtree_start = std::chrono::high_resolution_clock::now();
          rtrees_[i][j].range_search(Polygon<uint64_t>(Point<uint64_t>(begin_[0], begin_[1]), Point<uint64_t>(end_[0], end_[1])), answ);
          // auto rtree_end = std::chrono::high_resolution_clock::now();
          // rtree_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(rtree_end - rtree_start).count();
          // rtree_cnt += 1;
          total_answ.insert(total_answ.end(), answ.begin(), answ.end());
        }
#else
        total_cnt += by;
        total_cnt += ey;
#endif
      }
      total_cnt += total_answ.size();
      // for (int i = 0; i < total_answ.size(); i++)
      // {
      //   uint64_t x = total_answ[i]->get_Pmin().get_X();
      //   uint64_t y = total_answ[i]->get_Pmin().get_Y();
      //   if (x < begin_[0] || y < begin_[1] || x > end_[0] || y > end_[1])
      //   {
      //     printf("error\n");
      //     exit(-1);
      //   }
      // }
      return total_answ;
    }

    // std::vector<Polygon<uint64_t> *> range_searchX(KeyType begin, KeyType end)
    // {
    //   uint64_t begin_ = round(begin * ratio_);
    //   uint64_t end_ = round(end * ratio_);
    //   uint64_t min_y = round(min_keys_[1] * ratio_);
    //   uint64_t max_y = round(max_keys_[1] * ratio_);
    //   std::vector<Polygon<uint64_t> *> total_answ;
    //   uint_fast32_t bx = gf_.GetCoordX(begin_);
    //   uint_fast32_t ex = gf_.GetCoordX(end_);
    //   for (int i = bx; i <= ex; i++)
    //   {
    //     for (int j = 0; j < rtree_size_; j++)
    //     {
    //       std::vector<Polygon<uint64_t> *> answ;
    //       // auto rtree_start = std::chrono::high_resolution_clock::now();
    //       rtrees_[i][j].range_search(Polygon<uint64_t>(Point<uint64_t>(begin_, min_y), Point<uint64_t>(end_, max_y)), answ);
    //       // auto rtree_end = std::chrono::high_resolution_clock::now();
    //       // rtree_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(rtree_end - rtree_start).count();
    //       rtree_cnt += 1;
    //       total_answ.insert(total_answ.end(), answ.begin(), answ.end());
    //     }
    //   }
    //   total_cnt += total_answ.size();
    //   return total_answ;
    // }

    // std::vector<Polygon<uint64_t> *> range_searchY(KeyType begin, KeyType end)
    // {
    //   uint64_t begin_ = round(begin * ratio_);
    //   uint64_t end_ = round(end * ratio_);
    //   uint64_t min_x = round(min_keys_[0] * ratio_);
    //   uint64_t max_x = round(max_keys_[0] * ratio_);
    //   std::vector<Polygon<uint64_t> *> total_answ;
    //   for (int i = 0; i < rtree_size_; i++)
    //   {
    //     uint_fast32_t by = gf_.GetCoordY(i, begin_);
    //     uint_fast32_t ey = gf_.GetCoordY(i, end_);
    //     for (int j = by; j <= ey; j++)
    //     {
    //       if (j < 0 || j >= rtree_size_)
    //       {
    //         continue;
    //       }
    //       std::vector<Polygon<uint64_t> *> answ;
    //       // auto rtree_start = std::chrono::high_resolution_clock::now();
    //       rtrees_[i][j].range_search(Polygon<uint64_t>(Point<uint64_t>(min_x, begin_), Point<uint64_t>(max_x, end_)), answ);
    //       // auto rtree_end = std::chrono::high_resolution_clock::now();
    //       // rtree_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(rtree_end - rtree_start).count();
    //       rtree_cnt += 1;
    //       total_answ.insert(total_answ.end(), answ.begin(), answ.end());
    //     }
    //   }
    //   total_cnt += total_answ.size();
    //   return total_answ;
    // }

    std::vector<Polygon<uint64_t> *> knn(std::vector<KeyType> key, int k_num)
    {
      // float width = (max_keys_[0] - min_keys_[0]) * sqrt(k_num / 1.0 / key_num_);
      // float height = (max_keys_[1] - min_keys_[1]) * sqrt(k_num / 1.0 / key_num_);
      float width = sqrt(k_num / 1.0 / key_num_);
      float height = sqrt(k_num / 1.0 / key_num_);
      std::vector<KeyType> begin_(2);
      std::vector<KeyType> end_(2);
      std::vector<Polygon<uint64_t> *> answ;
      while (true)
      {
        knn_query_cnt++;
        // begin_[0] = max(float(0), key[0] - width / 2);
        // begin_[1] = max(float(0), key[1] - height / 2);
        // end_[0] = key[0] + width / 2;
        // end_[1] = key[1] + height / 2;
        begin_[0] = max(float(0), key[0] - width);
        begin_[1] = max(float(0), key[1] - height);
        end_[0] = key[0] + width;
        end_[1] = key[1] + height;
        std::vector<Polygon<uint64_t> *> cur_answ = range_search(begin_, end_);
        if (cur_answ.size() >= k_num)
        {
          break;
        }
        // width *= sqrt(2);
        // height *= sqrt(2);
        width *= 2;
        height *= 2;
      }
      return answ;
    }

  private:
    std::vector<std::vector<uint64_t>> data_;
    std::vector<float> min_keys_;
    std::vector<float> max_keys_;
    std::vector<Polygon<uint64_t>> pols_;
    std::vector<std::vector<RTree<uint64_t>>> rtrees_;
    // std::unordered_map<uint64_t, SpatialIndex::ISpatialIndex *> rtrees_;
    // std::unordered_map<uint64_t, int64_t> indexIdentifier_;
    gf::GridFile<uint64_t> gf_;
    const size_t num_radix_bits_;
    size_t group_size_;
    size_t rs_group_size_;
    size_t cell_size_;
    size_t exp_;
    size_t ratio_;
    size_t key_num_;
    size_t rtree_size_;
  };
} // namespace mdm