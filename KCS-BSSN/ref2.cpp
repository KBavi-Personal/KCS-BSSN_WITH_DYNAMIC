#include "Ref2.h"
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

#include "Dijkstra.h"
#include "BFS.h"
#include <chrono>
#include <queue>
#include "SampleData.h"

static inline uint64_t edge_key(int a, int b) {
	if (a > b) std::swap(a, b);
	return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

Ref2::Ref2(SampleData* _sd) : sd(_sd) {}
Ref2::~Ref2() {}
static std::vector<std::unordered_set<int>> BuildSNAdj(const SampleData* sd)
{
	const int n = sd->bn->no_of_users;
	std::vector<std::unordered_set<int>> adj(n);

	for (int i = 0; i < sd->no_remining_users; ++i) {
		int u = sd->remining_users[i];
		if (u < 0) continue;
		if (sd->bn->users_vec[u]->isPruned) continue;

		// out
		for (const auto& nb : sd->bn->users_vec[u]->outNeighbors) {
			int v = nb.first;
			if (v < 0 || v >= n) continue;
			if (sd->bn->users_vec[v]->isPruned) continue;
			if (u == v) continue;
			adj[u].insert(v);
			adj[v].insert(u);
		}
		// in
		for (const auto& nb : sd->bn->users_vec[u]->inNeighbors) {
			int v = nb.first;
			if (v < 0 || v >= n) continue;
			if (sd->bn->users_vec[v]->isPruned) continue;
			if (u == v) continue;
			adj[u].insert(v);
			adj[v].insert(u);
		}
	}
	return adj;
}

static int edge_support(const std::vector<std::unordered_set<int>>& adj, int u, int v)
{
	const auto& Nu = adj[u];
	const auto& Nv = adj[v];
	const auto& small = (Nu.size() < Nv.size()) ? Nu : Nv;
	const auto& large = (Nu.size() < Nv.size()) ? Nv : Nu;

	int sup = 0;
	for (int w : small) {
		if (w == u || w == v) continue;
		if (large.find(w) != large.end()) ++sup;
	}
	return sup;
}


void Ref2::remove_user_from_remaining_users(int index) {
	const int last = sd->no_remining_users - 1;
	if (index < 0 || index > last) return; 
	int removed_id = sd->remining_users[index];
	sd->bn->users_vec[removed_id]->isPruned = true;

	sd->remining_users[index] = sd->remining_users[last];
	sd->no_remining_users--;
}

void Ref2::remove_pruned_user() {
	for (int i = 0; i < sd->no_remining_users; ) {
		if (sd->remining_users[i] == -1) {
			sd->remining_users[i] = sd->remining_users[sd->no_remining_users - 1];
			sd->no_remining_users--;
		}
		else {
			++i;
		}
	}
}

void Ref2::remove_pruned_POI() {
	for (int i = 0; i < sd->no_remining_POI; ) {
		if (sd->remining_POI[i] == -1) {
			sd->remining_POI[i] = sd->remining_POI[sd->no_remining_POI - 1];
			sd->no_remining_POI--;
		}
		else {
			++i;
		}
	}
}

void Ref2::Keyword_Based_Ref2(const std::vector<int>& Q)
{
	std::unordered_set<int> qset;
	qset.reserve(Q.size() * 2);
	for (int k : Q) qset.insert(k);

	int i_poi = 0;
	for (auto* p : sd->poi_vec) {
		bool has_key = false;
		for (int kk = 0; kk < p->num_of_keys; ++kk) {
			if (qset.find(p->keys[kk]) != qset.end()) { has_key = true; break; }
		}
		if (has_key) sd->remining_POI[i_poi++] = p->id;
		else p->isPruned = true;
	}
	sd->no_remining_POI = i_poi;
	const int oldN = sd->no_remining_users;
	for (int i = 0; i < oldN; ++i) {
		int u_id = sd->remining_users[i];
		auto* u = sd->bn->users_vec[u_id];

		bool all_pruned = true;
		const auto& chk = u->checkin_locations;

		if (chk.empty()) all_pruned = true;

		for (const auto& pf : chk) {
			if (!sd->poi_vec[pf.first]->isPruned) { all_pruned = false; break; }
		}

		if (all_pruned) {
			u->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}

void Ref2::Omega_Based_Ref2(int omega)
{

	const int oldN = sd->no_remining_users;
	for (int i = 0; i < oldN; ++i) {
		const int u_id = sd->remining_users[i];
		int total_freq = 0;

		for (const auto& pf : sd->bn->users_vec[u_id]->checkin_locations) {
			if (sd->poi_vec[pf.first]->isPruned) continue;
			total_freq += pf.second;
		}

		if (total_freq < omega) {
			sd->bn->users_vec[u_id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}

void Ref2::Pi_Based_Ref2(float pi)
{

	std::unordered_map<int, std::pair<int, int>> stat;
	stat.reserve(sd->no_remining_POI * 2);

	for (int i = 0; i < sd->no_remining_POI; ++i)
		stat.emplace(sd->remining_POI[i], std::make_pair(0, 0));


	for (int i = 0; i < sd->no_remining_users; ++i) {
		const int u_id = sd->remining_users[i];
		for (const auto& pf : sd->bn->users_vec[u_id]->checkin_locations) {
			const int poi = pf.first;
			if (sd->poi_vec[poi]->isPruned) continue;

			auto it = stat.find(poi);
			if (it == stat.end()) continue;

			it->second.first += 1;       
			it->second.second += pf.second;
		}
	}

	for (auto& kv : stat) {
		const int poi = kv.first;
		const int fcount = kv.second.first;
		const int fsum = kv.second.second;
		const double avg = (fcount > 0) ? (double)fsum / (double)fcount : 0.0;
		if (avg < (double)pi) sd->poi_vec[poi]->isPruned = true;
	}
	int w = 0;
	for (int i = 0; i < sd->no_remining_POI; ++i) {
		int poi = sd->remining_POI[i];
		if (!sd->poi_vec[poi]->isPruned) sd->remining_POI[w++] = poi;
	}
	sd->no_remining_POI = w;
	for (int i = 0; i < sd->no_remining_users; ++i) {
		int u_id = sd->remining_users[i];
		bool all_pruned = true;
		for (const auto& pf : sd->bn->users_vec[u_id]->checkin_locations) {
			if (!sd->poi_vec[pf.first]->isPruned) { all_pruned = false; break; }
		}
		if (all_pruned) {
			sd->bn->users_vec[u_id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}

void Ref2::Social_Distance_Based_Ref2(int d, int q_index)
{

	BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, q_index, d, dist_tmp, true);

	for (int i = 0; i < sd->no_remining_users; ++i) {
		const int id = sd->remining_users[i];
		if (id < 0) continue;
		if (id >= (int)dist_tmp.size() || dist_tmp[id] == -1) {
			sd->bn->users_vec[id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}

#include <sstream>
#include "CommStru.h"
void Ref2::Structural_Cohesiveness_Ref2(int k)
{
	if (k <= 2) return;

	const int n = sd->bn->no_of_users;
	const int thr = k - 2;

	auto adj = BuildSNAdj(sd);

	std::unordered_map<uint64_t, int> sup;
	sup.reserve((size_t)sd->no_remining_users * 8);

	std::deque<std::pair<int, int>> q;


	for (int u = 0; u < n; ++u) {
		for (int v : adj[u]) {
			if (u >= v) continue;
			int s = edge_support(adj, u, v);
			sup.emplace(edge_key(u, v), s);
			if (s < thr) q.push_back({ u,v });
		}
	}


	auto remove_edge = [&](int u, int v) {
		adj[u].erase(v);
		adj[v].erase(u);
		};

	while (!q.empty()) {
		auto [u, v] = q.front();
		q.pop_front();

		if (u < 0 || u >= n || v < 0 || v >= n) continue;
		if (adj[u].find(v) == adj[u].end()) continue;

		uint64_t uvk = edge_key(u, v);
		auto it = sup.find(uvk);
		if (it == sup.end() || it->second >= thr) continue;

		std::vector<int> common;
		common.reserve(std::min(adj[u].size(), adj[v].size()));
		const auto& small = (adj[u].size() < adj[v].size()) ? adj[u] : adj[v];
		const auto& large = (adj[u].size() < adj[v].size()) ? adj[v] : adj[u];

		for (int w : small) {
			if (w == u || w == v) continue;
			if (large.find(w) != large.end()) common.push_back(w);
		}

		remove_edge(u, v);

		for (int w : common) {
			uint64_t uw = edge_key(u, w);
			uint64_t vw = edge_key(v, w);

			auto itu = sup.find(uw);
			if (itu != sup.end() && --itu->second < thr) q.push_back({ std::min(u,w), std::max(u,w) });

			auto itv = sup.find(vw);
			if (itv != sup.end() && --itv->second < thr) q.push_back({ std::min(v,w), std::max(v,w) });
		}
	}

	for (int i = 0; i < sd->no_remining_users; ++i) {
		int id = sd->remining_users[i];
		if (id >= 0 && adj[id].empty()) {
			sd->bn->users_vec[id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}

void Ref2::Spatial_Distance_Ref2_cache(float sigma, const std::set<int>& v_p_prime_set, int q_index)
{

	const int N = sd->total_n_of_intersection_points;
	if (N <= 0) return;

	if (v_p_prime_set.empty()) {
		for (int i = 0; i < sd->no_remining_users; ++i) {
			int u_id = sd->remining_users[i];
			if (u_id < 0) continue;
			if (u_id == q_index) continue;
			sd->bn->users_vec[u_id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
		remove_pruned_user();
		return;
	}


	const float INF = std::numeric_limits<float>::infinity();
	std::unordered_map<int, std::vector<float>> distCache;
	distCache.reserve(v_p_prime_set.size() * 2);
	auto computeSSSP = [&](int srcIP) -> const std::vector<float>&{
		auto it = distCache.find(srcIP);
		if (it != distCache.end()) return it->second;

		std::vector<float> dist(N, INF);
		if (srcIP < 0 || srcIP >= N) {
			auto [insIt, _] = distCache.emplace(srcIP, std::move(dist));
			return insIt->second;
		}

		dist[srcIP] = 0.0f;
		using Node = std::pair<float, int>; 
		std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
		pq.push({ 0.0f, srcIP });

		while (!pq.empty()) {
			auto [du, u] = pq.top();
			pq.pop();
			if (du != dist[u]) continue; 

			const auto& neighbors = sd->roadNetworkAdjMtrx[u]; 
			for (const auto& [v, w] : neighbors) {
				const float nd = du + (float)w;
				if (nd < dist[v]) {
					dist[v] = nd;
					pq.push({ nd, v });
				}
			}
		}

		auto [insIt, _] = distCache.emplace(srcIP, std::move(dist));
		return insIt->second;
		};


	for (int i = 0; i < sd->no_remining_users; ++i) {
		const int u_id = sd->remining_users[i];
		if (u_id < 0) continue;
		if (u_id == q_index) continue;

		const auto& chk = sd->bn->users_vec[u_id]->checkin_locations;
		if (chk.empty()) continue;

		const double denom = (double)chk.size();
		const double cutoff_sum = (double)sigma * denom;

		bool ok = false;
		for (int q_p_id : v_p_prime_set) {
			const int srcIP = sd->poi_vec[q_p_id]->my_closest_intesect_point;
			const auto& dist = computeSSSP(srcIP);

			double total = 0.0;
			bool fail_this_qp = false;

			for (const auto& u_p : chk) {
				const int poi_id = u_p.first;
				const int dstIP = sd->poi_vec[poi_id]->my_closest_intesect_point;

				const float d = (dstIP >= 0 && dstIP < N) ? dist[dstIP] : INF;

				if (!std::isfinite(d)) { fail_this_qp = true; break; }

				total += (double)d;

	
				if (total > cutoff_sum) { fail_this_qp = true; break; }
			}

			if (!fail_this_qp) {
				ok = true;
				break;
			}
		}

		if (!ok) {
			sd->bn->users_vec[u_id]->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}

	remove_pruned_user();

}

std::set<int> Ref2::Reset()
{
	std::set<int> remain;
	for (int i = 0; i < sd->bn->no_of_users; ++i) {
		sd->bn->users_vec[i]->isPruned = false;
		remain.insert(sd->bn->users_vec[i]->id);
	}
	for (int i = 0; i < sd->no_of_POI; ++i) {
		sd->poi_vec[i]->isPruned = false;
	}
	return remain;
}

bool Ref2::HasInflunceToAllRemaining(Dijkstra* dijkstra, float theta, int source_index, int q_index)
{
	const int u = sd->remining_users[source_index];
	std::vector<float> inf_u_to = dijkstra->MaxInfluenceFromSource_All_OutNeighbors(u, theta);
	for (int j = 0; j < sd->no_remining_users; ++j) {
		const int v = sd->remining_users[j];
		if (v == q_index || u == v) continue;
		if (inf_u_to[v] < theta) return false;
	}
	return true;
}
void Ref2::Influence_Based_Ref2_dijkstra(Dijkstra* dijkstra, float theta, int q_index)
{
	bool changed = true;

	while (changed) {
		changed = false;

		std::vector<float> inf_q_to = dijkstra->MaxInfluenceFromSource_All_OutNeighbors(q_index, theta);

		std::vector<char> prune_flag(sd->bn->no_of_users, 0);

		for (int i = 0; i < sd->no_remining_users; ++i) {
			int u_id = sd->remining_users[i];
			if (u_id == q_index) continue;

			if (inf_q_to[u_id] < theta) {
				prune_flag[u_id] = 1;
			}
		}


		for (int i = 0; i < sd->no_remining_users;) {
			int u_id = sd->remining_users[i];
			if (prune_flag[u_id]) {
				sd->bn->users_vec[u_id]->isPruned = true;
				remove_user_from_remaining_users(i);
				changed = true;
			}
			else {
				++i;
			}
		}
	}

}


std::string  Ref2::print_remain_candidates(std::string label) {
	std::ostringstream result;
	result << " " << label << "\n";
	std::set<int> remain_users_id_set;
	std::set<int> remain_poi_id_set;
	for (int i = 0; i < sd->no_remining_users; i++) {
		remain_users_id_set.insert(sd->remining_users[i]);
	}
	for (int i = 0; i < sd->no_remining_POI; i++) {
		remain_poi_id_set.insert(sd->remining_POI[i]);
	}
	result << "\n user(" + std::to_string(remain_users_id_set.size()) + "): ";
	for (auto i : remain_users_id_set) {
		result << std::to_string(i) + "#";
	}
	result << "\n poi(" + std::to_string(remain_poi_id_set.size()) + "): ";
	for (auto i : remain_poi_id_set) {
		result << std::to_string(i) + "#";
	}
	result << "\n";
	return result.str();
}





void Ref2::Run_Ref2(const std::vector<int>& Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index, Sextuple<bool, bool, bool, bool, bool, bool, bool> m) {

	Dijkstra dijkstra(sd);  

	if (m.keyword && !sd->bn->users_vec[q_index]->isPruned) {
		Keyword_Based_Ref2(Q);
	}

	int before;
	while (sd->no_remining_users > 1 && !sd->bn->users_vec[q_index]->isPruned) {

		before = sd->no_remining_users;
		if (m.pi) {
			Pi_Based_Ref2((float)pi);  

			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (m.omega) {
			Omega_Based_Ref2(omega);
			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (m.social) {
			Social_Distance_Based_Ref2(d, q_index);
			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (m.structrual) {
			Structural_Cohesiveness_Ref2(k); 
			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (m.spacial) {
			std::set<int> v_p_prime_set;
			for (const auto& pf : sd->bn->users_vec[q_index]->checkin_locations) {
				if (!sd->poi_vec[pf.first]->isPruned) v_p_prime_set.insert(pf.first);
			}
			Spatial_Distance_Ref2_cache(sigma, v_p_prime_set, q_index); 
			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (m.influnce) {
			Influence_Based_Ref2_dijkstra(&dijkstra, theta, q_index); 
			if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
				break;
		}
		if (before == sd->no_remining_users) { 
			break;
		}
	}

	if (sd->bn->users_vec[q_index]->isPruned) sd->no_remining_users = 0;
	return;
}




void Ref2::Keyword_Based_Ref2_dynamic(const std::vector<int>& Q)
{
	std::unordered_set<int> qset;
	qset.reserve(Q.size() * 2);
	for (int k : Q) qset.insert(k);
	std::vector<int> remaining_poi_ids;

	int i_poi = 0;
	for(int i = 0; i < sd->no_remining_POI; ++i) {
		int p_id = sd->remining_POI[i];
		auto* p = sd->poi_vec[p_id];
		bool has_key = false;
		for (int kk = 0; kk < p->num_of_keys; ++kk) {
			if (qset.find(p->keys[kk]) != qset.end()) { has_key = true; break; }
		}
		if (has_key) remaining_poi_ids.push_back(p->id);
		else p->isPruned = true;
	}
	sd->no_remining_POI = remaining_poi_ids.size();
	for(int i = 0; i < sd->no_remining_POI; ++i) {
		sd->remining_POI[i] = remaining_poi_ids[i];
	}

	// prune users whose all checkins are pruned
	const int oldN = sd->no_remining_users;
	for (int i = 0; i < oldN; ++i) {
		int u_id = sd->remining_users[i];
		auto* u = sd->bn->users_vec[u_id];

		bool all_pruned = true; // <-- FIX: set per user
		const auto& chk = u->checkin_locations;

		if (chk.empty()) all_pruned = true;

		for (const auto& pf : chk) {
			if (!sd->poi_vec[pf.first]->isPruned) { all_pruned = false; break; }
		}

		if (all_pruned) {
			u->isPruned = true;
			sd->remining_users[i] = -1;
		}
	}
	remove_pruned_user();
}





void Ref2::Run_Ref_dynamic(const std::vector<int>& Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index) {

	Dijkstra dijkstra(sd);

	if (!sd->bn->users_vec[q_index]->isPruned) {
		Keyword_Based_Ref2_dynamic(Q);
	}

	int before;
	while (sd->no_remining_users > 1 && !sd->bn->users_vec[q_index]->isPruned) {

		before = sd->no_remining_users;

		Pi_Based_Ref2((float)pi); 

		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		Omega_Based_Ref2(omega); 
		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		Social_Distance_Based_Ref2(d, q_index); 
		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		Structural_Cohesiveness_Ref2(k); 
		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		// V_p' : q's unpruned POIs
		std::set<int> v_p_prime_set;
		for (const auto& pf : sd->bn->users_vec[q_index]->checkin_locations) {
			if (!sd->poi_vec[pf.first]->isPruned) v_p_prime_set.insert(pf.first);
		}
		Spatial_Distance_Ref2_cache(sigma, v_p_prime_set, q_index); 
		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		Influence_Based_Ref2_dijkstra(&dijkstra, theta, q_index); 
		if (sd->no_remining_users < 1 || sd->bn->users_vec[q_index]->isPruned)
			break;

		if (before == sd->no_remining_users) { // no change, stop
			break;
		}
	}

	if (sd->bn->users_vec[q_index]->isPruned) sd->no_remining_users = 0;

	return;
}

