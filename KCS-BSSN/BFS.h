#pragma once
#include <vector>

class User;

class BFS {
public:
	static void ComputeDistancesFromSourceCutoff(const std::vector<User*>& users, int src, int max_d, std::vector<int>& dist_out, bool respect_pruned = false);
};