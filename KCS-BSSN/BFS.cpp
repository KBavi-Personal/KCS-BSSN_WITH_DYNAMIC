#include "BFS.h"
#include "User.h"
#include <queue>

void BFS::ComputeDistancesFromSourceCutoff(const std::vector<User*>& users, int src, int max_d, std::vector<int>& dist_out, bool respect_pruned)
{
	const int n = (int)users.size();
	dist_out.assign(n, -1);

	if (src < 0 || src >= n || users[src] == nullptr) return;
	if (respect_pruned && users[src]->isPruned) return;

	std::queue<int> q;
	dist_out[src] = 0;
	q.push(src);

	while (!q.empty()) {
		const int u = q.front();
		q.pop();

		if (respect_pruned && users[u]->isPruned) continue;

		const int du = dist_out[u];
		if (max_d >= 0 && du >= max_d) {
			// Cutoff: do not expand beyond maxDepth
			continue;
		}

		auto relax = [&](int v) {
			if (v < 0 || v >= n || users[v] == nullptr) return;
			if (respect_pruned && users[v]->isPruned) return;
			if (dist_out[v] != -1) return;
			dist_out[v] = du + 1;
			q.push(v);
			};

		for (const auto& e : users[u]->outNeighbors) {
			relax(e.first);
		}

	}
}