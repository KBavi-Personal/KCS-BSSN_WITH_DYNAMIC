#pragma once
#include <vector>
#include <unordered_set>
struct CommStru {
    int q_index = -1;
    std::unordered_set<int> user_ids;
    std::unordered_set<int> poi_ids;
};
