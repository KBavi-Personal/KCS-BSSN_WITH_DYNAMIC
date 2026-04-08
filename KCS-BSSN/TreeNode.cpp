#include "TreeNode.h"
#include <cfloat>
#include <string>
#include "SampleData.h"

TreeNode::TreeNode(int _level, int _P_s_size, int _P_r_size, Tree* _tree, bool isLeaf, int pivot_user_index, int _total_n_keys)
{
	id = -1; 
	depth = _level;
	P_r_size = _P_r_size;
	P_s_size = _P_s_size;
	my_tree = _tree;
	is_leaf = isLeaf;
	isroot = false;
	myUserIndexForPivot = pivot_user_index;
	ub_sub = 0;
	ub_w_in = 0;
	ub_w_out = 0;
	min_max_dist_P_s = new float[2 * P_s_size];
	for (int i = 0; i < P_s_size; i++) {
		min_max_dist_P_s[2 * i] = FLT_MAX;
		min_max_dist_P_s[2 * i + 1] = -FLT_MAX;
	}

	total_n_keys = _total_n_keys;
	Keys_ub_Fsum.resize(_total_n_keys, 0.0);
	Keys_ub_Fmax.resize(_total_n_keys, 0.0);
}
TreeNode::~TreeNode()
{
	if (min_max_dist_P_s != NULL && P_s_size > 0) {
		delete[] min_max_dist_P_s;
		min_max_dist_P_s = NULL;
		P_s_size = 0;
	}
}

float Get_lb_Avg_Dist_user_p(SampleData* sd, User* u, POI* p)
{
	if (u->checkin_locations.size() == 0)return 0.0f;
	float sum_dist = 0;
	float max_dist_r = 0;
	float dist_r = 0;
	for (auto u_p_freq : u->checkin_locations) {
		float max_dist_r = 0.0f;
		for (int i = 0; i < sd->P_r_size; i++)
		{
			dist_r = abs(sd->poi_vec[u_p_freq.first]->dist_P_r[i] - p->dist_P_r[i]);
			max_dist_r = dist_r > max_dist_r ? dist_r : max_dist_r;
		}
		sum_dist += max_dist_r;
	}

	return sum_dist / u->checkin_locations.size();
}

double Get_lb_Avg_Dist_two_users(SampleData* sd, User* u, User* v)
{
	if (v->checkin_locations.size() == 0)return 0.0;
	double lb_avg_dist = 0;
	for (auto v_p_freq : v->checkin_locations)
		lb_avg_dist += Get_lb_Avg_Dist_user_p(sd, u, sd->poi_vec[v_p_freq.first]);
	lb_avg_dist = lb_avg_dist / v->checkin_locations.size();

	return lb_avg_dist;
}