#pragma once
#include <vector>
#include <cstdint>

struct BatchVisit {
    int user_id = -1;
    int poi_id = -1;
    int day = 0;
};

struct BatchUpdateResult {
    uint64_t batch_id = 0;
    std::vector<int> updated_users;
};

struct UserVisit {
    int user_id = -1;
    int poi_id = -1;
    std::vector<int> days;
};



