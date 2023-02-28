#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <unordered_map>
#include "include/2D/multiD_map.h"
#include <stx/btree_multimap.h>
#include "libmorton/morton.h"
#include "include/kdtree/Tree.h"
using namespace std;
uint64_t rtree_duration = 0;
uint64_t spline_duration = 0;
uint64_t rtree_cnt = 0;
std::vector<Polygon<uint64_t>> pols_;
uint64_t x_cell_size;
uint64_t y_cell_size;
uint64_t minx;
uint64_t miny;
int k_num;
uint64_t t_size;
uint32_t max_prefix_all = 0;
unordered_map<uint64_t, uint64_t> cntid;
uint64_t avggrid_rtree_rangesearch(RTree<uint64_t> ***rtree, vector<float> data, vector<float> upper_data, size_t ratio)
{
  std::vector<uint64_t> plow(2);
  plow[0] = round(data[0] * ratio);
  plow[1] = round(data[1] * ratio);
  std::vector<uint64_t> phigh(2);
  phigh[0] = round(upper_data[0] * ratio);
  phigh[1] = round(upper_data[1] * ratio);
  std::vector<Polygon<uint64_t> *> total_answ;
  int bx, ex, by, ey;
  if (plow[0] < minx)
  {
    bx = 0;
  }
  else
  {
    bx = min(t_size - 1, (plow[0] - minx) / x_cell_size);
  }
  if (phigh[0] < minx)
  {
    ex = 0;
  }
  else
  {
    ex = min(t_size - 1, (phigh[0] - minx) / x_cell_size);
  }
  if (plow[1] < miny)
  {
    by = 0;
  }
  else
  {
    by = min(t_size - 1, (plow[1] - miny) / y_cell_size);
  }
  if (phigh[1] < miny)
  {
    ey = 0;
  }
  else
  {
    ey = min(t_size - 1, (phigh[1] - miny) / y_cell_size);
  }
  for (int i = bx; i <= ex; i++)
  {
    for (int j = by; j <= ey; j++)
    {
      std::vector<Polygon<uint64_t> *> answ;
      rtree[i][j]->range_search(Polygon<uint64_t>(Point<uint64_t>(plow[0], plow[1]), Point<uint64_t>(phigh[0], phigh[1])), answ);
      total_answ.insert(total_answ.end(), answ.begin(), answ.end());
    }
  }
  total_cnt += total_answ.size();
  // for (int i = 0; i < total_answ.size(); i++)
  // {
  //   uint64_t x = total_answ[i]->get_Pmin().get_X();
  //   uint64_t y = total_answ[i]->get_Pmin().get_Y();
  //   if (x < plow[0] || y < plow[1] || x > phigh[0] || y > phigh[1])
  //   {
  //     printf("error\n");
  //     exit(-1);
  //   }
  // }
  return total_answ.size();
}

uint64_t rstar_rangesearch(SpatialIndex::ISpatialIndex *tree, vector<float> data, vector<float> upper_data, size_t ratio)
{
  uint64_t xlow = round(data[0] * ratio), ylow = round(data[1] * ratio);
  uint64_t xhigh = round(upper_data[0] * ratio), yhigh = round(upper_data[1] * ratio);
  double plow[2], phigh[2];
  plow[0] = xlow;
  plow[1] = ylow;
  phigh[0] = xhigh;
  phigh[1] = yhigh;
  util::MyVisitor vis;
  SpatialIndex::Region r = SpatialIndex::Region(plow, phigh, 2);
  tree->intersectsWithQuery(r, vis);
  total_cnt += rstar_res.size();
  // for (uint64_t res : rstar_res)
  // {
  //   if (res < xlow || res > xhigh)
  //   {
  //     printf("error\n");
  //     // exit(-1);
  //   }
  // }
  uint64_t z = rstar_res.size();
  rstar_res.clear();
  return z;
}

uint64_t avggrid_btree_rangesearch(stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>, stx::btree_default_map_traits<uint64_t, uint64_t>> ***bt, vector<float> data, vector<float> upper_data, size_t ratio)
{
  std::vector<uint64_t> plow(2);
  plow[0] = round(data[0] * ratio);
  plow[1] = round(data[1] * ratio);
  std::vector<uint64_t> phigh(2);
  phigh[0] = round(upper_data[0] * ratio);
  phigh[1] = round(upper_data[1] * ratio);
  std::vector<Polygon<uint64_t> *> total_answ;
  int bx, ex, by, ey;
  if (plow[0] < minx)
  {
    bx = 0;
  }
  else
  {
    bx = min(t_size - 1, (plow[0] - minx) / x_cell_size);
  }
  if (phigh[0] < minx)
  {
    ex = 0;
  }
  else
  {
    ex = min(t_size - 1, (phigh[0] - minx) / x_cell_size);
  }
  if (plow[1] < miny)
  {
    by = 0;
  }
  else
  {
    by = min(t_size - 1, (plow[1] - miny) / y_cell_size);
  }
  if (phigh[1] < miny)
  {
    ey = 0;
  }
  else
  {
    ey = min(t_size - 1, (phigh[1] - miny) / y_cell_size);
  }
  int cur_cnt = 0;
  for (uint64_t i = bx; i <= ex; i++)
  {
    for (uint64_t j = by; j <= ey; j++)
    {
      std::vector<Polygon<uint64_t> *> answ;
      uint64_t x = max(i * x_cell_size + minx, plow[0]);
      uint64_t y = max(j * y_cell_size + miny, plow[1]);
      uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
      auto iter = bt[i][j]->lower_bound(id);
      uint64_t resx, resy;
      while (iter != bt[i][j]->end())
      {
        id = iter->first;
        libmorton::morton2D_64_decode(id, resx, resy);
        resx += i * x_cell_size + minx;
        resy += j * y_cell_size + miny;
        if (resx >= plow[0] && resx <= phigh[0] && resy >= plow[1] && resy <= phigh[1])
        {
          total_cnt++;
          cur_cnt++;
        }
        iter++;
      }
    }
  }
  return cur_cnt;
}

int main(int argc, char **argv)
{
  if (argc != 13)
  {
    cerr << "usage: " << argv[0] << "<data_file> <data_num> <query_num> <dimension> <exp> <index_type> <num_radix_bits> <cell_size> <rs_group_size> <query_type> (<range_size>/<k_num>) (<update_ratio>)" << endl;
    throw;
  }
  srand(921224);
  const string data_file = argv[1];
  const int data_size = atoi(argv[2]);
  const int query_size = atoi(argv[3]);
  const int dim = atoi(argv[4]);
  const int exp = atoi(argv[5]);
  string index_type = argv[6];
  const int num_radix_bits = atoi(argv[7]);
  const int cell_size = atoi(argv[8]);
  const int rs_group_size = atoi(argv[9]);
  string query_type = argv[10];
  int range_size;
  double update_ratio;
  if (query_type == "range")
  {
    range_size = atoi(argv[11]);
  }
  if (query_type == "update")
  {
    range_size = atoi(argv[11]);
    update_ratio = atof(argv[12]);
  }
  if (query_type == "knn")
  {
    k_num = atoi(argv[11]);
  }
  // const int cluster_num = atoi(argv[13]);
  size_t ratio = pow(10, exp);
  // vector<vector<float>> srcdata = util::LoadData<float>(data_file, dim, data_size);
  vector<vector<float>> data = util::LoadData<float>(data_file, dim);
  std::random_shuffle(data.begin(), data.end());
  data.resize(data_size);
  int n = data.size();
  // int dima_num = pow(n / cell_size, 0.75) * cell_size;
  // int dimb_num = pow(n / cell_size, 0.5) * cell_size;
  // std::sort(srcdata.begin(), srcdata.end(), [](const std::vector<float> &a, const std::vector<float> &b)
  //           { return a[0] < b[0]; });
  // std::sort(srcdata.begin(), srcdata.begin() + dima_num, [](const std::vector<float> &a, const std::vector<float> &b)
  //           { return a[1] < b[1]; });

  // vector<vector<float>> data;
  // for (int i = 0; i < dima_num; i++)
  // {
  //   vector<float> item(2);
  //   item[0] = srcdata[i][2];
  //   item[1] = srcdata[i][3];
  //   data.push_back(item);
  // }

  t_size = std::sqrt(n / cell_size) + 1;
  vector<float> min_keys = util::GetMinKey<float>(data, dim);
  vector<float> max_keys = util::GetMaxKey<float>(data, dim);
  util::gridsize[0] = (max_keys[0] - min_keys[0]) / range_size;
  util::gridsize[1] = (max_keys[1] - min_keys[1]) / range_size;
  util::gridsizeX = (max_keys[0] - min_keys[0]) / range_size / range_size;
  util::gridsizeY = (max_keys[1] - min_keys[1]) / range_size / range_size;
  x_cell_size = uint64_t((max_keys[0] - min_keys[0]) * ratio) / (t_size - 2);
  y_cell_size = uint64_t((max_keys[1] - min_keys[1]) * ratio) / (t_size - 2);
  minx = round(min_keys[0] * ratio);
  miny = round(min_keys[1] * ratio);

  vector<vector<float>> insert_data;
  if (query_type == "update")
  {
    insert_data.insert(insert_data.end(), data.begin() + int(update_ratio * data_size), data.end());
    data.resize(int(update_ratio * data_size));
    cout << data.size() << endl;
  }
  std::sort(data.begin(), data.end(), [](const std::vector<float> &a, const std::vector<float> &b)
            { return a[0] < b[0]; });
  n = data.size();
  int key_num = n;

  // grid+r-tree
  mdm::MultiDMap<float> *map;

  // kd-tree
  Tree *kd_tree = new Tree();

  // avggrid+r-tree
  RTree<uint64_t> ***rtree = (RTree<uint64_t> ***)malloc(sizeof(RTree<uint64_t> **) * t_size);
  for (int i = 0; i < t_size; i++)
  {
    rtree[i] = (RTree<uint64_t> **)malloc(sizeof(RTree<uint64_t> *) * t_size);
  }

  // r*-tree
  SpatialIndex::IStorageManager *memoryfile = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  SpatialIndex::StorageManager::IBuffer *file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*memoryfile, 10, false);
  int64_t rtree_id;
  SpatialIndex::ISpatialIndex *tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, rtree_id);

  // avggrid+b-tree
  stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                      stx::btree_default_map_traits<uint64_t, uint64_t>>
      ***bt = (stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                   stx::btree_default_map_traits<uint64_t, uint64_t>> ***)malloc(sizeof(stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                                                                                                            stx::btree_default_map_traits<uint64_t, uint64_t>> **) *
                                                                                                 t_size);
  for (int i = 0; i < t_size; i++)
  {
    bt[i] = (stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                 stx::btree_default_map_traits<uint64_t, uint64_t>> **)malloc(sizeof(stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                                                                                                         stx::btree_default_map_traits<uint64_t, uint64_t>> *) *
                                                                                              t_size);
  }

  auto build_begin = chrono::high_resolution_clock::now();

  if (index_type == "grid+r-tree")
  {
    // grid+r-tree
    map = new mdm::MultiDMap<float>(data, num_radix_bits, cell_size, rs_group_size, exp);
  }

  if (index_type == "avggrid+r-tree")
  {
    // avggrid+r-tree
    for (int i = 0; i < t_size; i++)
    {
      for (int j = 0; j < t_size; j++)
      {
        rtree[i][j] = new RTree<uint64_t>(100);
      }
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(x, y)));
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = pols_[i].get_Pmin().get_X(), y = pols_[i].get_Pmin().get_Y();
      rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert_polygon(&pols_[i], &pols_[i]);
    }
  }

  if (index_type == "r*-tree")
  {
    // r*-tree
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      std::string value = std::to_string(x);
      double p[2];
      p[0] = x;
      p[1] = y;
      SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
      tree->insertData((uint32_t)(value.size() + 1), reinterpret_cast<const uint8_t *>(value.c_str()), r, x);
    }
  }

  if (index_type == "avggrid+b-tree")
  {
    // avggrid+b-tree
    for (int i = 0; i < t_size; i++)
    {
      for (int j = 0; j < t_size; j++)
      {
        bt[i][j] = new stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                           stx::btree_default_map_traits<uint64_t, uint64_t>>();
      }
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
      if (cntid.find(id) != cntid.end())
      {
        cntid[id] += 1;
      }
      else
      {
        cntid[id] = 1;
      }
      bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert2(id, x);
      // cnt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]++;
    }
  }

  if (index_type == "kd-tree")
  {
    // kd-tree
    for (int i = 0; i < n; i++)
    {
      kd_tree->insertObject(data[i], i + 1);
    }
  }

  auto build_end = chrono::high_resolution_clock::now();
  double build_s = chrono::duration_cast<chrono::nanoseconds>(build_end - build_begin).count() / 1.0 / 1e9;
  cout << "build/s: " << build_s << endl;

  std::random_shuffle(data.begin(), data.end());
  data.resize(query_size);
  vector<vector<float>> upper_data = util::GetUpperData<float>(data, dim, max_keys);
  vector<float> upper_dataX = util::GetUpperDataX<float>(data, dim, max_keys);
  vector<float> upper_dataY = util::GetUpperDataY<float>(data, dim, max_keys);
  n = data.size();
  int testcnt = 3;
  uint_fast64_t k = 0;
  // knn_query_cnt = 0;
  spline_duration = 0;
  auto lookup_begin = chrono::high_resolution_clock::now();
  for (int test = 0; test < testcnt; test++)
  {
    if (query_type == "range" || query_type == "update")
    {
      // range query

      if (index_type == "grid+r-tree")
      {
        // grid+r-tree

        for (int i = 0; i < n; i++)
        {
          map->range_search(data[i], upper_data[i]);
        }

        // for (int i = 0; i < n; i++)
        // {
        //   map->range_searchX(data[i][0], upper_dataX[i]);
        // }

        // for (int i = 0; i < n; i++)
        // {
        //   map->range_searchY(data[i][1], upper_dataY[i]);
        // }
      }

      if (index_type == "avggrid+r-tree")
      {
        // avggrid+r-tree
        for (int i = 0; i < n; i++)
        {
          k |= avggrid_rtree_rangesearch(rtree, data[i], upper_data[i], ratio);
        }
      }

      if (index_type == "r*-tree")
      {
        // r*-tree
        for (int i = 0; i < n; i++)
        {
          k |= rstar_rangesearch(tree, data[i], upper_data[i], ratio);
        }
      }

      if (index_type == "avggrid+b-tree")
      {
        // avggrid+b-tree
        for (int i = 0; i < n; i++)
        {
          avggrid_btree_rangesearch(bt, data[i], upper_data[i], ratio);
        }
      }

      if (index_type == "kd-tree")
      {
        // kd-tree
        for (int i = 0; i < n; i++)
        {
          std::vector<uint32_t> results = kd_tree->rangeSearch(data[i], upper_data[i]);
          // printf("%d\n",results.size());
          total_cnt += results.size();
        }
      }
    }
    else if (query_type == "knn")
    {
      // knn query

      if (index_type == "grid+r-tree")
      {
        // grid+r-tree
        for (int i = 0; i < n; i++)
        {
          map->knn(data[i], k_num);
        }
      }

      if (index_type == "avggrid+r-tree")
      {
        // avggrid+r-tree
        for (int i = 0; i < n; i++)
        {
          float width = sqrt(k_num / 1.0 / key_num);
          float height = sqrt(k_num / 1.0 / key_num);
          std::vector<float> begin_(2);
          std::vector<float> end_(2);
          while (true)
          {
            knn_query_cnt++;
            begin_[0] = max(float(0), data[i][0] - width);
            begin_[1] = max(float(0), data[i][1] - height);
            end_[0] = data[i][0] + width;
            end_[1] = data[i][1] + height;
            uint64_t answ_size = avggrid_rtree_rangesearch(rtree, begin_, end_, ratio);
            if (answ_size >= k_num)
            {
              break;
            }
            width *= 2;
            height *= 2;
          }
        }
      }

      if (index_type == "r*-tree")
      {
        // r*-tree
        for (int i = 0; i < n; i++)
        {
          float width = sqrt(k_num / 1.0 / key_num);
          float height = sqrt(k_num / 1.0 / key_num);
          std::vector<float> begin_(2);
          std::vector<float> end_(2);
          while (true)
          {
            knn_query_cnt++;
            begin_[0] = max(float(0), data[i][0] - width);
            begin_[1] = max(float(0), data[i][1] - height);
            end_[0] = data[i][0] + width;
            end_[1] = data[i][1] + height;
            uint64_t answ_size = rstar_rangesearch(tree, begin_, end_, ratio);
            if (answ_size >= k_num)
            {
              break;
            }
            width *= 2;
            height *= 2;
          }
        }
      }

      if (index_type == "avggrid+b-tree")
      {
        // avggrid+b-tree
        for (int i = 0; i < n; i++)
        {
          float width = sqrt(k_num / 1.0 / key_num);
          float height = sqrt(k_num / 1.0 / key_num);
          std::vector<float> begin_(2);
          std::vector<float> end_(2);
          while (true)
          {
            knn_query_cnt++;
            begin_[0] = max(float(0), data[i][0] - width);
            begin_[1] = max(float(0), data[i][1] - height);
            end_[0] = data[i][0] + width;
            end_[1] = data[i][1] + height;
            uint64_t answ_size = avggrid_btree_rangesearch(bt, begin_, end_, ratio);
            if (answ_size >= k_num)
            {
              break;
            }
            width *= 2;
            height *= 2;
          }
        }
      }
    }
    else
    {
      // point query

      if (index_type == "grid+r-tree")
      {
        // grid+r-tree
        for (int i = 0; i < n; i++)
        {
          k |= map->find(data[i]);
        }
      }

      if (index_type == "avggrid+r-tree")
      {
        // avggrid+r-tree
        for (int i = 0; i < n; i++)
        {
          uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
          std::vector<Polygon<uint64_t> *> answ;
          rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->range_search(Polygon<uint64_t>(Point<uint64_t>(x, y)), answ);
          total_cnt += answ.size();
          // if (answ[0]->get_Pmin().get_X() != x || answ[0]->get_Pmin().get_Y() != y)
          // {
          //   printf("error result\n");
          //   exit(-1);
          // }
          k |= answ[0]->get_Pmin().get_X();
        }
      }

      if (index_type == "r*-tree")
      {
        // r*-tree
        for (int i = 0; i < n; i++)
        {
          uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
          double p[2];
          p[0] = x;
          p[1] = y;
          util::MyVisitor vis;
          SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
          tree->intersectsWithQuery(r, vis);
          total_cnt += rstar_res.size();
          // if (rstar_res[0] != x)
          // {
          //   printf("error\n");
          //   // exit(-1);
          // }
          k |= rstar_res[0];
          rstar_res.clear();
        }
      }

      if (index_type == "avggrid+b-tree")
      {
        // avggrid+b-tree
        for (int i = 0; i < n; i++)
        {
          uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
          uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
          auto iter = bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->lower_bound(id);
          uint64_t res = iter->second;
          // if (res != x)
          // {
          //   printf("error\n");
          //   // exit(-1);
          // }
          while (iter->first == id)
          {
            total_cnt++;
            iter++;
          }
          k |= res;
        }
      }
    }
  }
  auto lookup_end = chrono::high_resolution_clock::now();
  uint64_t lookup_ns = chrono::duration_cast<chrono::nanoseconds>(lookup_end - lookup_begin).count();
  cout << "testcnt: " << testcnt << endl;
  cout << "lookup/ns: " << lookup_ns / n / testcnt << endl;
  cout << "total_cnt: " << total_cnt / testcnt << endl;
  cout << "knn_query_cnt/bs_distance: " << knn_query_cnt / testcnt << endl;
  cout << "spline_duration: " << spline_duration / testcnt << endl;
  cout << k << endl;

  if (query_type == "update")
  {
    auto insert_begin = chrono::high_resolution_clock::now();
    n = insert_data.size();
    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      // for (int i = 0; i < n; i++)
      // {
      //   map->insert(insert_data[i]);
      // }
    }
    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(x, y)));
      }
      for (int i = 0; i < n; i++)
      {
        uint64_t x = pols_[i].get_Pmin().get_X(), y = pols_[i].get_Pmin().get_Y();
        rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert_polygon(&pols_[i], &pols_[i]);
      }
    }
    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        std::string value = std::to_string(x);
        double p[2];
        p[0] = x;
        p[1] = y;
        SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
        tree->insertData((uint32_t)(value.size() + 1), reinterpret_cast<const uint8_t *>(value.c_str()), r, x);
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
        bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert2(id, x);
        // cnt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]++;
      }
    }

    auto insert_end = chrono::high_resolution_clock::now();
    uint64_t insert_ns = chrono::duration_cast<chrono::nanoseconds>(insert_end - insert_begin).count();
    cout << "insert/ns: " << insert_ns / n << endl;

    n = data.size();
    total_cnt = 0;
    auto afterlookup_begin = chrono::high_resolution_clock::now();
    k = 0;

    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      for (int i = 0; i < n; i++)
      {
        map->range_search(data[i], upper_data[i]);
      }
    }

    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        k |= avggrid_rtree_rangesearch(rtree, data[i], upper_data[i], ratio);
      }
    }

    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        k |= rstar_rangesearch(tree, data[i], upper_data[i], ratio);
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        avggrid_btree_rangesearch(bt, data[i], upper_data[i], ratio);
      }
    }

    auto afterlookup_end = chrono::high_resolution_clock::now();
    uint64_t afterlookup_ns = chrono::duration_cast<chrono::nanoseconds>(afterlookup_end - afterlookup_begin).count();
    cout << "afterlookup/ns: " << afterlookup_ns / n << endl;
    cout << "total_cnt: " << total_cnt << endl;
    cout << k << endl;
  }
  // cout << "rtree_duration: " << rtree_duration / testcnt << endl;
  // cout << "rtree_cnt: " << rtree_cnt / testcnt << endl;
  // cout << "per_rtree: " << rtree_duration / rtree_cnt << endl;

  return 1;
}