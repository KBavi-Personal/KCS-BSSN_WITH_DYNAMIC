#include "Pruning.h"
#include <cmath>
#include <cfloat>  

Pruning::Pruning(SampleData* _sd)
{
	sd = _sd;
}
Pruning::~Pruning()
{
}

bool Pruning::Keyword_based_pruning(User* u, std::vector<int>& Q) const 
{
	for (int key : Q)
		if (u->Keys_visited_count[key] > 0)
			return false;
	return true;
}

bool Pruning::Omega_based_pruning(User* u, int Omega)  const
{
	if (u->ub_f_sum < Omega) {
		return true;
	}
	return false;
}

bool Pruning::Pi_based_pruning(User* u, float pi) const
{
	if (u->ub_f_avg < pi) {
		return true;
	}
	return false;
}


bool Pruning::Influnce_based_pruning(User* q, User* v, float theta) const
{
	if (q == v) return false;
	float ub_ISF_q_v, ub_ISF_v_q;
	if (q->outNeighbors.count(v->id) == 1) 
	{
		ub_ISF_q_v = fmax(q->ub_w_out, v->ub_w_in);
	
	}
	else
	{
		ub_ISF_q_v = q->ub_w_out * v->ub_w_in;
		
	}

	if (ub_ISF_q_v < theta){
		return true;
	}
	return false;
}


bool Pruning::Structural_cohesiveness_pruning_sub(User* u, int k) const
{
	if (u->ub_sub < (k - 2)) {
		return true;
	}
	return false;
}

bool Pruning::Social_distance_based_pruning(User* u, User* q, int d) const
{
	if (!u || !q) return true;
	if (d < 0) return true;

	int lb_dist_s = 0;
	bool has_valid_pivot = false;

	for (int i = 0; i < sd->P_s_size; ++i) {
		const int du = u->dist_P_s[i];
		const int dq = q->dist_P_s[i];

		
		if (du < 0 || dq < 0) continue;

		has_valid_pivot = true;
		const int diff = std::abs(du - dq);
		if (diff > lb_dist_s) lb_dist_s = diff;

	
		if (lb_dist_s > d)
			return true;
	}

	
	if (!has_valid_pivot)
		return false;

	return false;
}

float Pruning::LB_poi_poi(int poiA, int poiB) const
{
	float lb = 0.0f;
	for (int i = 0; i < sd->P_r_size; ++i) {
		float d = std::fabs(sd->poi_vec[poiA]->dist_P_r[i] - sd->poi_vec[poiB]->dist_P_r[i]);
		if (d > lb) lb = d;
	}
	return lb;
}
float Pruning::Get_lb_avgdist(User* u, User* q, std::vector<int>& Q) const
{
	float lb_avgdist = FLT_MAX;
	float lb_avgdist_u_p = 0;
	float sum_dist = 0;
	float max_dist_r = 0;
	float dist_r = 0;

	if (u->checkin_locations.size() == 0) return lb_avgdist;
	for (auto q_p_freq : q->checkin_locations) {
		auto p = sd->poi_vec[q_p_freq.first];
		if (p->isPruned) continue;
		sum_dist = 0;
		if (p->HasAnyKey(Q)) {
			for (auto u_p_freq : u->checkin_locations) {
				sum_dist += LB_poi_poi(u_p_freq.first, q_p_freq.first);
			}
			lb_avgdist_u_p = sum_dist / (float)u->checkin_locations.size();

			lb_avgdist = lb_avgdist_u_p < lb_avgdist ? lb_avgdist_u_p : lb_avgdist;
		}
	}
	return lb_avgdist;
}

bool Pruning::Spatial_distance_based_pruning(User* u, User* q, std::vector<int>& Q, float sigma) const
{
	float lb_avgdist = Get_lb_avgdist(u, q, Q);
	if (lb_avgdist > sigma) {
		return true;
	}
	return false;
}

bool Pruning::Keyword_based_pruning_index_node(TreeNode* node, std::vector<int>& Q, int d) const
{
	for (int query_key : Q)

		if (node->Keys_ub_Fsum[query_key] > 0) {
			return false;
		}
	return true;
}

bool Pruning::Omega_based_pruning_index_node(TreeNode* node, std::vector<int>& Q, int Omega, int d) const
{
	double max_ub_f_sum = 0;
	for (int query_key : Q)
	{
		max_ub_f_sum = std::max(node->Keys_ub_Fsum[query_key], max_ub_f_sum);
	}
	if (max_ub_f_sum < Omega) {
		return true;
	}
	return false;
}

bool Pruning::Pi_based_pruning_index_node(TreeNode* node, std::vector<int>& Q, float Pi, int d) const
{
	double max_ub_f_max = 0;
	double ub_f_max = 0;
	for (int query_key : Q)
	{
		max_ub_f_max = std::max(node->Keys_ub_Fmax[query_key], max_ub_f_max);
	}
	if (max_ub_f_max < Pi) {
		return true;
	}
	return false;
}

bool Pruning::Influnce_based_pruning_index_node(TreeNode* node, User* q, float theta) const
{
	float ub_ISF_q_node = q->ub_w_out * node->ub_w_in;
	
	if (ub_ISF_q_node < theta) {
		return true;
	}
	return false;
}

bool Pruning::Structural_cohesiveness_pruning_sub_index_node(TreeNode* node, int k) const
{
	if (node->ub_sub < (k - 2)) {
		return true;
	}
	return false;
}


bool Pruning::Social_distance_based_pruning_index_node(TreeNode* node, User* q, int d) const
{

	int lb_dist_s = 0;
	for (int i = 0; i < sd->P_s_size; i++) {
		int dq = q->dist_P_s[i];
		int dmin = node->min_max_dist_P_s[2 * i];
		int dmax = node->min_max_dist_P_s[2 * i + 1];

		int dist_to_interval = 0;
		if (dq < dmin) dist_to_interval = dmin - dq;
		else if (dq > dmax) dist_to_interval = dq - dmax;

		lb_dist_s = std::max(lb_dist_s, dist_to_interval);
		if (lb_dist_s > d) return true;
	}
	return false;
}