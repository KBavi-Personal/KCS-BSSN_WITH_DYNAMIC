#include "OfflineCalculations.h"
#include "Dijkstra.h"
#include "MinHeap.h"
#include "Triple.h"
#include <iostream>
#include"BFS.h"
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <cfloat>
#include <cmath>

OfflineCalculations::OfflineCalculations(SampleData* _sd)
{
	sd = _sd;
}

OfflineCalculations::~OfflineCalculations()
{
}
#include <sstream>
void OfflineCalculations::Run(int pivot_iterations, std::string filename, int distribution_type)
{
	Calculate_edges_sup_in_SN2();
	sd->larget_f_sum = 0;
	sd->larget_f_avg = 0;
	Compute_POI_fsum_count_per_key();
	Compute_poi_max_freq();

	for (User* u : sd->bn->users_vec) {
		//Calculate_u_keys(u);
		u->ub_sub = Get_ub_sub(u);
		u->ub_w_in = Get_ub_w_in(u);
		u->ub_w_out = Get_ub_w_out(u);
		u->ub_f_sum = Get_ub_f_sum(u);
		sd->larget_f_sum = sd->larget_f_sum < u->ub_f_sum ? u->ub_f_sum : sd->larget_f_sum;
		u->ub_f_avg = Get_ub_f_avg(u);
	}
	sd->larget_f_avg = Get_largest_f_avg();
	Compute_Key_fsum_fmax_for_users();
	SN_Pivot_Selection(pivot_iterations);
	RN_Pivot_Selection(pivot_iterations);
}
int OfflineCalculations::Get_ub_f_sum(User* u)
{
	int ub_f_sum = 0;
	for (std::pair<int, int> p : u->checkin_locations)
		ub_f_sum += p.second;
	return ub_f_sum;
}
void OfflineCalculations::Compute_poi_max_freq()
{
	const int P_N = sd->no_of_POI;
	poi_max_freq.assign(P_N, 0);

	
	for (int uid = 0; uid < (int)sd->bn->users_vec.size(); ++uid) {
		User* u = sd->bn->users_vec[uid];
		if (!u) continue;

		for (const auto& pf : u->checkin_locations) {
			const int poi = pf.first;
			const int f = pf.second;
			if (poi < 0 || poi >= P_N) continue;

			if (f > poi_max_freq[poi]) poi_max_freq[poi] = f;
		}
	}
}
float OfflineCalculations::Get_ub_f_avg(User* u)
{
	if (!u || u->checkin_locations.empty()) return 0.0f;

	int ub = 0;
	for (const auto& pf : u->checkin_locations) {
		const int poi = pf.first;
		if (poi < 0 || poi >= (int)poi_max_freq.size()) continue;
		ub = std::max(ub, poi_max_freq[poi]);
	}
	return (float)ub;

}

float OfflineCalculations::Get_largest_f_avg()
{
	std::map<int, std::pair<int, int>> poi_id_fcount_fsum;
	float largest_f_avg = 0;
	for (User* u : sd->bn->users_vec)
		for (auto p_freq : u->checkin_locations) {
			poi_id_fcount_fsum[p_freq.first].first += 1;
			poi_id_fcount_fsum[p_freq.first].second += p_freq.second;
		}

	for (auto p : poi_id_fcount_fsum) {
		float avg_freq = 0;
		if (p.second.first > 0)
			avg_freq = p.second.second / p.second.first;
		largest_f_avg = avg_freq > largest_f_avg ? avg_freq : largest_f_avg;
	}
	return largest_f_avg;
}
float OfflineCalculations::Get_ub_w_in(User* u)
{
	float ub_w_in = 0.0;
	for (std::pair<int, float> inNeighbor : u->inNeighbors)
		if (inNeighbor.second > ub_w_in)
			ub_w_in = inNeighbor.second;
	return ub_w_in;
}
float OfflineCalculations::Get_ub_w_out(User* u)
{
	float ub_w_out = 0.0;
	for (std::pair<int, float> outNeighbor : u->outNeighbors)
		if (outNeighbor.second > ub_w_out)
			ub_w_out = outNeighbor.second;
	return ub_w_out;
}
float OfflineCalculations::Get_ub_sub(User* u)
{
	float ub_sub_u = 0.0;
	std::pair<int, int> e;
	
	for (std::pair<int, float> v : u->inNeighbors) {
		float sub_e = (float)sd->bn->Get_edge_sup(v.first, u->id);
		if (sub_e > ub_sub_u)
			ub_sub_u = sub_e;
	}
	for (std::pair<int, float> v : u->outNeighbors) {
		float sub_e = (float)sd->bn->Get_edge_sup(v.first, u->id);
		if (sub_e > ub_sub_u)
			ub_sub_u = sub_e;
	}
	return ub_sub_u;
}

bool IsUserInVec(const std::vector<User*>& userVec, int id) {
	for (int i = 0; i < userVec.size(); i++)
	{
		if (userVec[i]->id == id)return true;
	}
	return false;
}


bool IsInVec(std::vector<int> vec, int index) {
	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] == index)return true;
	}
	return false;
}


void OfflineCalculations::Print_SN_Pivot(std::string filename)
{

	std::string results_to_print;
	for (int j = 0; j < sd->P_s_size; j++) {

		results_to_print += std::to_string(sd->P_s_vec[j]->id) + " ";
	}
	results_to_print += "\n";
	for (int i = 0; i < sd->bn->no_of_users; i++)
	{
		results_to_print += std::to_string(i) + " ";
		for (int j = 0; j < sd->P_s_size; j++) {

			results_to_print += std::to_string(sd->bn->users_vec[i]->dist_P_s[j]) + " ";
		}
		results_to_print += "\n";
	}
	sd->Print_To_File(filename + ".txt", results_to_print);
}
void OfflineCalculations::Print_RN_Pivot(std::string filename)
{

	std::string results_to_print;
	for (int j = 0; j < sd->P_r_size; j++) {

		results_to_print += std::to_string(sd->P_r_vec[j]) + " ";
	}
	results_to_print += "\n";
	for (int i = 0; i < sd->no_of_POI; i++)
	{
		results_to_print += std::to_string(i) + " ";
		for (int j = 0; j < sd->P_r_size; j++) {

			results_to_print += std::to_string(sd->poi_vec[i]->dist_P_r[j]) + " ";
		}
		results_to_print += "\n";
	}
	sd->Print_To_File(filename + ".txt", results_to_print);
}

static int intersection_size_sorted(const std::vector<int>& a, const std::vector<int>& b)
{
	int i = 0, j = 0, cnt = 0;
	while (i < (int)a.size() && j < (int)b.size()) {
		if (a[i] == b[j]) { ++cnt; ++i; ++j; }
		else if (a[i] < b[j]) ++i;
		else ++j;
	}
	return cnt;
}

void OfflineCalculations::Calculate_edges_sup_in_SN2()
{
	sd->bn->Reset_edges_sup_values();

	const int n = sd->bn->no_of_users;
	std::vector<std::vector<int>> adj(n);

	for (const auto& e : sd->bn->Edges_sup()) {
		int u = e.first.first;
		int v = e.first.second;
		if (u == v) continue;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	for (int i = 0; i < n; ++i) {
		auto& vec = adj[i];
		std::sort(vec.begin(), vec.end());
		vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
	}

	for (auto& e : sd->bn->Edges_sup_mut()) {
		const int u = e.first.first;
		const int v = e.first.second;
		e.second = intersection_size_sorted(adj[u], adj[v]);
	}

}


static std::vector<std::vector<int>>
BFS_users_by_distance(const std::vector<User*>& user_vec, int src, int d_max)
{
	const int n = (int)user_vec.size();
	std::vector<int> dist(n, -1);
	std::queue<int> q;

	dist[src] = 0;
	q.push(src);

	std::vector<std::vector<int>> buckets(d_max + 1);

	while (!q.empty()) {
		const int u = q.front();
		q.pop();

		const int du = dist[u];
		if (du == d_max) continue;

	
		auto relax = [&](int v) {
			if (v < 0 || v >= n) return;
			if (dist[v] != -1) return;

			dist[v] = du + 1;
			if (dist[v] <= d_max) {
				buckets[dist[v]].push_back(v);
				q.push(v);
			}
			};

		for (const auto& kv : user_vec[u]->outNeighbors) {
			relax(kv.first);
		}
	
	}
	return buckets;
}


void OfflineCalculations::Compute_Key_fsum_fmax_for_users()
{
	
	for (int i = 0; i < sd->bn->no_of_users; ++i) {
		User* user = sd->bn->users_vec[i];

		for (const auto& p_pair : user->checkin_locations) {
			POI* p = sd->poi_vec[p_pair.first];
			const int c = p_pair.second;

		
			for (int kk = 0; kk < p->num_of_keys; ++kk) {
				const int key = p->keys[kk];
				user->Keys_Fsum[key] += c; 
				user->Keys_Fmax[key] = std::max(user->Keys_Fmax[key], (double)c);; 
				user->Keys_visited_count[key] += 1; 
			}
		}
	}
}

static inline bool finite_dist(float x) {
	return x != FLT_MAX && std::isfinite(x);
}


static std::vector<double> Compute_POI_cost_sums(SampleData* sd)
{
	const int P = sd->no_of_POI;
	const int R = sd->P_r_size;

	std::vector<double> S(P, 0.0);

	std::vector<std::vector<float>> D(P);
	D.reserve(P);
	for (int i = 0; i < P; ++i) {
		D[i] = sd->poi_vec[i]->dist_P_r;
	}

	for (int a = 0; a < P; ++a) {
		const std::vector<float> Da = D[a];
		double sum = 0.0;

		for (int p = 0; p < P; ++p) {
			if (p == a) continue;

			const std::vector<float> Dp = D[p];

			float mx = 0.0f;
			bool any = false;

			for (int j = 0; j < R; ++j) {
				const float xa = Da[j];
				const float xp = Dp[j];
				if (!finite_dist(xa) || !finite_dist(xp)) continue; 
				const float diff = std::fabs(xa - xp);
				if (diff > mx) mx = diff;
				any = true;
			}

			if (any) sum += mx;
			
		}

		S[a] = sum;
	}

	return S;
}

static double Compute_Total_cost_from_S(SampleData* sd, const std::vector<double>& S)
{
	double total = 0.0;

	for (User* u : sd->bn->users_vec) {
		for (const auto& p_freq : u->checkin_locations) {
			const int c = p_freq.first;

			total += S[c];
		}
	}
	return total;
}

void OfflineCalculations::RN_Pivot_Selection(int iteration)
{
	sd->P_r_vec.clear();
	sd->P_r_vec.reserve(sd->P_r_size);

	Dijkstra dijkstra(sd);


	std::unordered_set<int> chosen;
	chosen.reserve(sd->P_r_size * 2);

	while ((int)sd->P_r_vec.size() < sd->P_r_size) {
		int ip = random1::uniform_int_pos(0, sd->total_n_of_intersection_points - 1);
		if (chosen.insert(ip).second) {
			sd->P_r_vec.push_back(ip);
		
			dijkstra.Calculate_shotest_dist_road_s_to_all((int)sd->P_r_vec.size() - 1);
		}
	}

	
	std::vector<double> S = Compute_POI_cost_sums(sd);
	double best_cost = Compute_Total_cost_from_S(sd, S);

	
	std::vector<float> old_pivot_dist_to_poi_vec;
	old_pivot_dist_to_poi_vec.reserve(sd->no_of_POI);

	for (int it = 0; it < iteration; ++it) {

	
		if (it == iteration - 1) break;

	
		const int old_pivot_index = random1::uniform_int_pos(0, (int)sd->P_r_vec.size() - 1);
		const int old_ip = sd->P_r_vec[old_pivot_index];

	
		int new_ip = random1::uniform_int_pos(0, sd->total_n_of_intersection_points - 1);
		while (chosen.count(new_ip)) {
			new_ip = random1::uniform_int_pos(0, sd->total_n_of_intersection_points - 1);
		}

		
		old_pivot_dist_to_poi_vec.clear();
		for (POI* p : sd->poi_vec) {
			old_pivot_dist_to_poi_vec.push_back(p->dist_P_r[old_pivot_index]);
		}

	
		chosen.erase(old_ip);
		chosen.insert(new_ip);
		sd->P_r_vec[old_pivot_index] = new_ip;

		
		dijkstra.Calculate_shotest_dist_road_s_to_all(old_pivot_index);

		
		std::vector<double> S_new = Compute_POI_cost_sums(sd);
		double new_cost = Compute_Total_cost_from_S(sd, S_new);

		if (new_cost < best_cost) {
		
			best_cost = new_cost;
			S.swap(S_new);
		}
		else {
			
			chosen.erase(new_ip);
			chosen.insert(old_ip);
			sd->P_r_vec[old_pivot_index] = old_ip;

			for (int i = 0; i < sd->no_of_POI; ++i) {
				sd->poi_vec[i]->dist_P_r[old_pivot_index] = old_pivot_dist_to_poi_vec[i];
			}
		}
	}
}


static int pick_reachable_pivot_uid(SampleData* sd, const std::vector<char>& isPivot, std::vector<int>& dist_tmp) {
	const int U = sd->bn->no_of_users;
	int best_uid = -1;
	int best_reach = -1;

	for (int t = 0; t < 500; ++t) {
		int uid = random1::uniform_int_pos(0, U - 1);
		while (isPivot[uid]) uid = random1::uniform_int_pos(0, U - 1);

		BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, uid, -1, dist_tmp, false);

		int reach = 0;
		for (int i = 0; i < U; ++i) {
			if (dist_tmp[i] >= 0) ++reach; 
		}

		if (reach > best_reach) {
			best_reach = reach;
			best_uid = uid;
		}
	}

	return best_uid;
}

static unsigned long long sum_pairwise_abs(std::vector<int>& vals)
{
	std::sort(vals.begin(), vals.end());
	unsigned long long total = 0;
	unsigned long long prefix = 0;

	for (size_t i = 0; i < vals.size(); ++i) {
		total += (unsigned long long)vals[i] * i - prefix;
		prefix += vals[i];
	}
	return total;
}

void OfflineCalculations::SN_Pivot_Selection(int iteration)
{
	const int U = sd->bn->no_of_users;
	const int P = sd->P_s_size;

	sd->P_s_vec.clear();
	sd->P_s_vec.reserve(P);
	std::vector<int> dist_tmp;



	std::vector<char> isPivot(U, 0);

	
	for (int j = 0; j < P; ++j) {
		int uid = pick_reachable_pivot_uid(sd, isPivot, dist_tmp);

		isPivot[uid] = 1;
		sd->P_s_vec.push_back(sd->bn->users_vec[uid]);


		BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, sd->P_s_vec[j]->id, -1, dist_tmp, false);
		for (int i = 0; i < sd->bn->no_of_users; ++i) {
			sd->bn->users_vec[i]->dist_P_s[j] = dist_tmp[i];
		}

	}

	auto compute_cost_proxy = [&]() -> unsigned long long {
		unsigned long long cost = 0;

		
		std::vector<int> vals;
		vals.reserve(U);

		for (int j = 0; j < P; ++j) {
			vals.clear();

			const int pivotId = sd->P_s_vec[j]->id;
			for (int i = 0; i < U; ++i) {
				const int uid = sd->bn->users_vec[i]->id;
				if (uid == pivotId) continue;

				const int d = sd->bn->users_vec[i]->dist_P_s[j];
				if (d == -1) continue; 
				vals.push_back(d);
			}

			if (vals.size() >= 2) cost += sum_pairwise_abs(vals);
		}
		return cost;
		};

	unsigned long long best_cost = compute_cost_proxy();


	std::vector<int> old_col;
	old_col.reserve(U);

	for (int it = 0; it < iteration; ++it) {
		if (it == iteration - 1) break;

		
		const int old_pivot_index = random1::uniform_int_pos(0, P - 1);
		const int old_uid = sd->P_s_vec[old_pivot_index]->id;

		
		int new_uid = pick_reachable_pivot_uid(sd, isPivot, dist_tmp);


	
		old_col.clear();
		for (int i = 0; i < U; ++i) old_col.push_back(sd->bn->users_vec[i]->dist_P_s[old_pivot_index]);

		
		isPivot[old_uid] = 0;
		isPivot[new_uid] = 1;
		sd->P_s_vec[old_pivot_index] = sd->bn->users_vec[new_uid];

	
		BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, sd->P_s_vec[old_pivot_index]->id, -1, dist_tmp, false);
		for (int i = 0; i < sd->bn->no_of_users; ++i) {
			sd->bn->users_vec[i]->dist_P_s[old_pivot_index] = dist_tmp[i];
		}

		const unsigned long long new_cost = compute_cost_proxy();

		if (new_cost < best_cost) {
			best_cost = new_cost; 
		}
		else {
		
			isPivot[new_uid] = 0;
			isPivot[old_uid] = 1;
			sd->P_s_vec[old_pivot_index] = sd->bn->users_vec[old_uid];

			for (int i = 0; i < U; ++i) sd->bn->users_vec[i]->dist_P_s[old_pivot_index] = old_col[i];
		}
	}
}

void OfflineCalculations::Compute_POI_fsum_count_per_key()
{
	for (User* u : sd->bn->users_vec)
	{
		if (!u) continue;
		for (const auto& pf : u->checkin_locations)
		{
			int poi_id = pf.first;
			POI* p = sd->poi_vec[poi_id];

			for (int i = 0; i < p->num_of_keys; i++) {
				p->Keys_Fsum[p->keys[i]] += pf.second;
				p->Keys_visited_count[p->keys[i]] += 1;
			}

		}

	}
}
