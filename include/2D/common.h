#pragma once

#include <cstddef>
#include <cstdint>
#include <spatialindex/SpatialIndex.h>
#include <immintrin.h>
std::vector<uint64_t> rstar_res;
int total_cnt=0;
int knn_query_cnt=0;
namespace rs
{
    template <class KeyType>
    struct Coord
    {
        KeyType x;
        double y;
    };

    struct FitParam
    {
        double slope;
        double intercept;
    };

} // namespace rs
