#include "Indexing.h"
#include <iostream>
#include "Dijkstra.h"
#include <cmath>
#include "rand.h"
#include "Sextuple.h"
#include <cfloat>
#include <vector>
#include <limits>
#include <algorithm>

Indexing::Indexing(SampleData* _sd)
{
	sd = _sd;
}
Indexing::~Indexing()
{
}
float Indexing::Get_Avg_Dist(User* u, POI* p)
{
	if (u->checkin_locations.size() == 0)return 0.0f;
	float sum_dist = 0;
	float max_dist_r = 0;
	float dist_r = 0;
	for (auto u_p_freq : u->checkin_locations) {
		max_dist_r = 0;
		for (int i = 0; i < sd->P_r_size; i++)
		{
			dist_r = abs(sd->poi_vec[u_p_freq.first]->dist_P_r[i] - p->dist_P_r[i]);
			max_dist_r = dist_r > max_dist_r ? dist_r : max_dist_r;
		}
		sum_dist += max_dist_r;
	}

	return sum_dist / u->checkin_locations.size();
}
std::pair<double, double> Indexing::Get_bs_score(User* u, User* v)
{
	double fsum = 0.0, favg = 0.0;
	int count = 0;
	if (u->checkin_locations.size() == 0 || v->checkin_locations.size() == 0) 	return { fsum, favg };;
	for (int i = 0; i < sd->toatl_n_of_used_keys; ++i) {
		count = 0;
		if (u->Keys_visited_count[i] > 0 && v->Keys_visited_count[i] > 0) {
			for (auto p_f : u->checkin_locations) {
				POI* p = sd->poi_vec[p_f.first];
				fsum += p->Keys_Fsum[i];
				favg += p->Keys_Fsum[i] / p->Keys_visited_count[i];
				count++;
			}
			if (count > 0) {
				favg = favg / count;
				fsum = fsum / count;
			}
		}
	}
	return { fsum, favg };
}

double Indexing::Get_rs_score(User* u, User* v)
{
	if (v->checkin_locations.size() == 0)return 0.0;
	double rs_score = 0;
	for (auto v_p_freq : v->checkin_locations)
		rs_score += Get_Avg_Dist(u, sd->poi_vec[v_p_freq.first]);
	rs_score = rs_score / v->checkin_locations.size();

	//printf("rs_score= %f \t", rs_score);
	return rs_score;
}
Triple<double, double, double> Indexing::Get_ss_score(User* u, User* v)
{
	Triple<double, double, double> ss_score{ 0,0,0 };
	ss_score.first = u->ub_sub + v->ub_sub;

	double ub_ISF;
	if (u->outNeighbors.count(v->id) == 1) //if e(u,v) in Edges
		ub_ISF = fmax(u->ub_w_out, v->ub_w_in);
	else
		ub_ISF = u->ub_w_out * v->ub_w_in;

	ss_score.second = ub_ISF;
	double lb_dist_s = 0;
	double dist = 0;
	for (int i = 0; i < sd->P_s_size; i++)
	{
		dist = abs(u->dist_P_s[i] - v->dist_P_s[i]);
		lb_dist_s = dist > lb_dist_s ? dist : lb_dist_s;
	}
	ss_score.third = lb_dist_s;
	return ss_score;
}
std::pair<double, double> Indexing::Get_bs_score_node(TreeNode* N, User* v)
{
	int fsum = 0, fmax = 0;
	for (int i = 0; i < sd->toatl_n_of_used_keys; i++) {
		if (N->Keys_ub_Fsum[i] > 0 && v->Keys_visited_count[i] > 0) {
			fsum += N->Keys_ub_Fsum[i];
			fsum += v->Keys_Fsum[i];
			fmax += N->Keys_ub_Fmax[i];
			if (v->Keys_visited_count[i] > 0)
				fmax += v->Keys_Fsum[i] / v->Keys_visited_count[i];
		}
	}

	return std::make_pair<double, double>(fsum, fmax);
}
Triple<double, double, double> Indexing::Get_ss_score_node(TreeNode* N, User* v)
{
	Triple<double, double, double> ss_score{ 0,0,0 };
	ss_score.first = N->ub_sub + v->ub_sub;

	float dist_min = 0, dist_max = 0;
	float lb_dist_s = 0, lb_dist_min = 0, lb_dist_max = 0;
	for (int i = 0; i < sd->P_s_size; i++)
	{

		dist_min = abs(v->dist_P_s[i] - N->min_max_dist_P_s[2 * i]);
		lb_dist_min = dist_min > lb_dist_min ? dist_min : lb_dist_min;

		dist_max = abs(v->dist_P_s[i] - N->min_max_dist_P_s[2 * i + 1]);
		lb_dist_max = dist_max > lb_dist_max ? dist_max : lb_dist_max;
	}
	lb_dist_s = lb_dist_min < lb_dist_max ? lb_dist_min : lb_dist_max;


	float ub_ISF;
	if (lb_dist_s <= 1) { //probably there is an edge between v and a user in N
		ub_ISF = fmax(N->ub_w_out, v->ub_w_in);
		lb_dist_s = 1;
	}
	else
		ub_ISF = N->ub_w_out * v->ub_w_in;
	ss_score.second = ub_ISF;
	ss_score.third = lb_dist_s;
	return ss_score;
}
void Indexing::Indexing_Pivot_Selection(int iteration, int piv_size, float bs_w, float ss_w, float rs_w)
{
	piv_size = std::max(1, piv_size);
	int user_index;
	float p_index_cost = FLT_MAX, temp_p_index_cost = 0;
	std::set<User*> temp_P_index_set;

	for (int i = 0; i < piv_size; i++)
	{
		user_index = random1::uniform_int_pos(0, sd->bn->no_of_users - 1);
		while (temp_P_index_set.count(sd->bn->users_vec[user_index]) == 1)//check if it already in list
			user_index = random1::uniform_int_pos(0, sd->bn->no_of_users - 1);
		temp_P_index_set.insert(sd->bn->users_vec[user_index]);
	}
	P_index_vec.clear();
	for (User* piv : temp_P_index_set)
		P_index_vec.push_back(piv);
	Partition_SN(bs_w, ss_w, rs_w);
	bool doCost = true;
	for (int i = 0; i < iteration; i++)
	{
		if (doCost)
			p_index_cost = calculate_p_index_cost(bs_w, ss_w, rs_w);
		doCost = false;
		//select a new pivot
		user_index = random1::uniform_int_pos(0, sd->bn->no_of_users - 1);
		while (temp_P_index_set.count(sd->bn->users_vec[user_index]) == 1)
			user_index = random1::uniform_int_pos(0, sd->bn->no_of_users - 1);

		int temp_index = random1::uniform_int_pos(0, P_index_vec.size() - 1);
		int temp_user_index = P_index_vec[temp_index]->id;

		temp_P_index_set.erase(P_index_vec[temp_index]);
		temp_P_index_set.insert(sd->bn->users_vec[user_index]);



		P_index_vec.clear();
		for (User* piv : temp_P_index_set)
			P_index_vec.push_back(piv);
		std::map < User*, std::vector<User*>> temp_Subgraphs = Subgraphs;// pivot_index, subgraph
		Partition_SN(bs_w, ss_w, rs_w);
		temp_p_index_cost = calculate_p_index_cost(bs_w, ss_w, rs_w);
		if (temp_p_index_cost > p_index_cost)
		{//back to orignal version of pivots
			temp_P_index_set.erase(sd->bn->users_vec[user_index]);
			temp_P_index_set.insert(sd->bn->users_vec[temp_user_index]);
			P_index_vec.clear();
			for (User* piv : temp_P_index_set)
				P_index_vec.push_back(piv);
			Subgraphs = temp_Subgraphs;
		}
		else
			p_index_cost = temp_p_index_cost;

	}
}
void Indexing::Partition_SN(float bs_w, float ss_w, float rs_w)
{
	Subgraphs.clear();
	if (P_index_vec.empty()) return;

	auto norm01 = [](double x, double mn, double mx) -> double {
		if (mx <= mn) return 0.0;
		double v = (x - mn) / (mx - mn);
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		return v;
		};

	struct Raw6 {
		double rs, bs_sum, bs_favg, ss_sum, ss_ub, ss_lb;
	};

	std::vector<Raw6> raw;
	raw.reserve(P_index_vec.size());
	const double EPS = 1e-12; // tie tolerance
	for (User* u : sd->bn->users_vec)
	{
		raw.clear();

		double mn_rs = DBL_MAX, mn_bs_sum = DBL_MAX, mn_bs_favg = DBL_MAX, mn_ss_sum = DBL_MAX, mn_ss_ub = DBL_MAX, mn_ss_lb = DBL_MAX;
		double mx_rs = 0.0, mx_bs_sum = 0.0, mx_bs_favg = 0.0, mx_ss_sum = 0.0, mx_ss_ub = 0.0, mx_ss_lb = 0.0;

		for (User* piv : P_index_vec)
		{
			Raw6 r{};
			r.rs = Get_rs_score(u, piv);

			auto fsum_favg = Get_bs_score(u, piv);
			r.bs_sum = fsum_favg.first;
			r.bs_favg = fsum_favg.second;

			auto ss = Get_ss_score(u, piv);
			r.ss_sum = ss.first;
			r.ss_ub = ss.second;
			r.ss_lb = ss.third;

			raw.push_back(r);

			mx_rs = std::max(mx_rs, r.rs);
			mx_bs_sum = std::max(mx_bs_sum, r.bs_sum);
			mx_bs_favg = std::max(mx_bs_favg, r.bs_favg);
			mx_ss_sum = std::max(mx_ss_sum, r.ss_sum);
			mx_ss_ub = std::max(mx_ss_ub, r.ss_ub);
			mx_ss_lb = std::max(mx_ss_lb, r.ss_lb);

			mn_rs = std::min(mn_rs, r.rs);
			mn_bs_sum = std::min(mn_bs_sum, r.bs_sum);
			mn_bs_favg = std::min(mn_bs_favg, r.bs_favg);
			mn_ss_sum = std::min(mn_ss_sum, r.ss_sum);
			mn_ss_ub = std::min(mn_ss_ub, r.ss_ub);
			mn_ss_lb = std::min(mn_ss_lb, r.ss_lb);
		}

		User* best_pivot = nullptr;
		double best_quality = -DBL_MAX;

		for (size_t i = 0; i < P_index_vec.size(); ++i)
		{
			const Raw6& r = raw[i];

			const double n_rs = norm01(r.rs, mn_rs, mx_rs);

			const double n_bs_sum = norm01(r.bs_sum, mn_bs_sum, mx_bs_sum);
			const double n_bs_favg = norm01(r.bs_favg, mn_bs_favg, mx_bs_favg);
			const double bs = 0.5 * (n_bs_sum + n_bs_favg);

			const double n_ss_sum = norm01(r.ss_sum, mn_ss_sum, mx_ss_sum);
			const double n_ss_ub = norm01(r.ss_ub, mn_ss_ub, mx_ss_ub);
			const double n_ss_lb = norm01(r.ss_lb, mn_ss_lb, mx_ss_lb);

			const double ss = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;

			double rs_sim = 1.0 - n_rs;
			double quality = (bs_w * bs) + (ss_w * ss) + (rs_w * rs_sim); // higher is better

			User* piv = P_index_vec[i];
			if (quality > best_quality) {
				best_quality = quality;
				best_pivot = piv;
			}
			else if (std::abs(quality - best_quality) < EPS) {
				best_quality = quality + EPS;
				best_pivot = piv;
			}
		}
		if (!best_pivot) best_pivot = P_index_vec.front();
		Subgraphs[best_pivot].push_back(u);

	}
}

double Indexing::calculate_p_index_cost(float bs_w, float ss_w, float rs_w)
{
	if (P_index_vec.empty() || sd == nullptr || sd->bn == nullptr) return 0.0;
	if (sd->bn->users_vec.empty()) return 0.0;

	auto norm01 = [](double x, double mn, double mx) -> double {
		if (mx <= mn) return 0.0;
		double v = (x - mn) / (mx - mn);
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		return v;
		};

	struct Raw6 { double rs, bs_sum, bs_max, ss_sum, ss_ub, ss_lb; };

	std::vector<Raw6> raw;
	raw.reserve(P_index_vec.size());

	double total_cost = 0.0;

	for (User* u : sd->bn->users_vec)
	{
		raw.clear();

		double mn_rs = DBL_MAX, mn_bs_sum = DBL_MAX, mn_bs_max = DBL_MAX, mn_ss_sum = DBL_MAX, mn_ss_ub = DBL_MAX, mn_ss_lb = DBL_MAX;
		double mx_rs = 0.0, mx_bs_sum = 0.0, mx_bs_max = 0.0, mx_ss_sum = 0.0, mx_ss_ub = 0.0, mx_ss_lb = 0.0;

		for (User* piv : P_index_vec)
		{
			Raw6 r{};
			r.rs = Get_rs_score(u, piv);

			auto fsum_fmax = Get_bs_score(u, piv);
			r.bs_sum = fsum_fmax.first;
			r.bs_max = fsum_fmax.second;

			auto ss = Get_ss_score(u, piv);
			r.ss_sum = ss.first;
			r.ss_ub = ss.second;
			r.ss_lb = ss.third;

			raw.push_back(r);

			mx_rs = std::max(mx_rs, r.rs);
			mx_bs_sum = std::max(mx_bs_sum, r.bs_sum);
			mx_bs_max = std::max(mx_bs_max, r.bs_max);
			mx_ss_sum = std::max(mx_ss_sum, r.ss_sum);
			mx_ss_ub = std::max(mx_ss_ub, r.ss_ub);
			mx_ss_lb = std::max(mx_ss_lb, r.ss_lb);

			mn_rs = std::min(mn_rs, r.rs);
			mn_bs_sum = std::min(mn_bs_sum, r.bs_sum);
			mn_bs_max = std::min(mn_bs_max, r.bs_max);
			mn_ss_sum = std::min(mn_ss_sum, r.ss_sum);
			mn_ss_ub = std::min(mn_ss_ub, r.ss_ub);
			mn_ss_lb = std::min(mn_ss_lb, r.ss_lb);
		}

		double best = -DBL_MAX;
		std::vector<double> score(P_index_vec.size(), 0.0);

		for (size_t i = 0; i < P_index_vec.size(); ++i)
		{
			const Raw6& r = raw[i];

			const double n_rs = norm01(r.rs, mn_rs, mx_rs);

			const double n_bs_sum = norm01(r.bs_sum, mn_bs_sum, mx_bs_sum);
			const double n_bs_max = norm01(r.bs_max, mn_bs_max, mx_bs_max);
			const double bs = 0.5 * (n_bs_sum + n_bs_max);

			const double n_ss_sum = norm01(r.ss_sum, mn_ss_sum, mx_ss_sum);
			const double n_ss_ub = norm01(r.ss_ub, mn_ss_ub, mx_ss_ub);
			const double n_ss_lb = norm01(r.ss_lb, mn_ss_lb, mx_ss_lb);
			const double ss = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;

			const double sc = (bs_w * bs) + (ss_w * ss) + (rs_w * n_rs);
			score[i] = sc;
			best = std::max(best, sc);
		}

		for (double sc : score) total_cost += (best - sc);
	}

	return total_cost;
}


std::vector<User*> Indexing::Indexing_Pivot_Selection_node(int iteration, const std::vector<TreeNode*>& Nodes, const std::vector<User*>& P_index_vec, int piv_size, float bs_w, float ss_w)
{
	const int M = (int)P_index_vec.size();
	if (M == 0 || piv_size <= 0) return {};
	if (piv_size > M) piv_size = M;

	std::vector<int> chosen;
	chosen.reserve(piv_size);
	std::vector<char> inSet(M, 0);

	while ((int)chosen.size() < piv_size) {
		const int idx = random1::uniform_int_pos(0, M - 1);
		if (!inSet[idx]) {
			inSet[idx] = 1;
			chosen.push_back(idx);
		}
	}

	auto buildPivotVector = [&]() {
		std::vector<User*> pivs;
		pivs.reserve(chosen.size());
		for (int idx : chosen) pivs.push_back(P_index_vec[idx]);
		return pivs;
		};

	float best_cost = calculate_p_index_cost_node(Nodes, buildPivotVector(), bs_w, ss_w);

	if (iteration <= 0 || piv_size == M) {
		return buildPivotVector();
	}

	for (int it = 0; it < iteration; ++it) {
		const int pos = random1::uniform_int_pos(0, (int)chosen.size() - 1);
		const int current_idx = chosen[pos];

		int cand_idx = random1::uniform_int_pos(0, M - 1);
		while (inSet[cand_idx]) cand_idx = random1::uniform_int_pos(0, M - 1);

		// swap
		inSet[current_idx] = 0;
		inSet[cand_idx] = 1;
		chosen[pos] = cand_idx;

		const float new_cost = calculate_p_index_cost_node(Nodes, buildPivotVector(), bs_w, ss_w);

		if (new_cost <= best_cost) {
			best_cost = new_cost;
		}
		else {
			// rollback
			inSet[cand_idx] = 0;
			inSet[current_idx] = 1;
			chosen[pos] = current_idx;
		}
	}

	return buildPivotVector();
}

float Indexing::calculate_p_index_cost_node(const std::vector<TreeNode*>& Nodes, const std::vector<User*>& pivots, float bs_w, float ss_w)
{
	if (Nodes.empty() || pivots.empty()) return 0.0f;

	struct Raw {
		double bs_sum; // fsum
		double bs_max; // fmax
		double ss_sum; // sum_sub
		double ss_ub;  // ub_ISF
		double ss_lb;  // lb_dist_s
	};

	auto norm01 = [](double x, double mn, double mx) -> double {
		if (mx <= mn) return 0.0;
		double v = (x - mn) / (mx - mn);
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		return v;
		};

	std::vector<Raw> raw(pivots.size());
	std::vector<double> score(pivots.size());

	double total_cost = 0.0;

	for (TreeNode* N : Nodes)
	{
		Raw mx{ 0,0,0,0,0 };
		Raw mn{ DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX };

		// 1) raw + min/max
		for (size_t i = 0; i < pivots.size(); ++i)
		{
			User* piv = pivots[i];

			const auto fsum_fmax = Get_bs_score_node(N, piv);
			const auto ss_score = Get_ss_score_node(N, piv);

			Raw r{
				(double)fsum_fmax.first,
				(double)fsum_fmax.second,
				(double)ss_score.first,
				(double)ss_score.second,
				(double)ss_score.third
			};

			raw[i] = r;

			mx.bs_sum = std::max(mx.bs_sum, r.bs_sum);
			mx.bs_max = std::max(mx.bs_max, r.bs_max);
			mx.ss_sum = std::max(mx.ss_sum, r.ss_sum);
			mx.ss_ub = std::max(mx.ss_ub, r.ss_ub);
			mx.ss_lb = std::max(mx.ss_lb, r.ss_lb);

			mn.bs_sum = std::min(mn.bs_sum, r.bs_sum);
			mn.bs_max = std::min(mn.bs_max, r.bs_max);
			mn.ss_sum = std::min(mn.ss_sum, r.ss_sum);
			mn.ss_ub = std::min(mn.ss_ub, r.ss_ub);
			mn.ss_lb = std::min(mn.ss_lb, r.ss_lb);
		}

		double best = -DBL_MAX;

		for (size_t i = 0; i < pivots.size(); ++i)
		{
			const Raw& r = raw[i];

			const double n_bs_sum = norm01(r.bs_sum, mn.bs_sum, mx.bs_sum);
			const double n_bs_max = norm01(r.bs_max, mn.bs_max, mx.bs_max);
			const double bs = 0.5 * (n_bs_sum + n_bs_max);

			const double n_ss_sum = norm01(r.ss_sum, mn.ss_sum, mx.ss_sum);
			const double n_ss_ub = norm01(r.ss_ub, mn.ss_ub, mx.ss_ub);
			const double n_ss_lb = norm01(r.ss_lb, mn.ss_lb, mx.ss_lb);

			const double ss = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;

			const double sc = (bs_w * bs) + (ss_w * ss); // higher is better
			score[i] = sc;
			best = std::max(best, sc);
		}

		for (double sc : score) total_cost += (best - sc);
	}

	return (float)total_cost;
}

