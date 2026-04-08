#include "DynamicApproch.h"
#include "Tree.h"
#include "ref2.h"
#include "rand.h"
#include "TreeTextIO.h"
#include <chrono>
#include <sstream>
#include "OfflineCalculations.h"
#include <iostream>

DynamicApproch::DynamicApproch()
    : sd(nullptr), inx(nullptr)
{
}
DynamicApproch::~DynamicApproch()
{
}

void DynamicApproch::RecomputeLeafMetrics(TreeNode* leaf)
{
	if (!leaf) return;
	std::fill(leaf->Keys_ub_Fsum.begin(), leaf->Keys_ub_Fsum.end(), 0.0);
	std::fill(leaf->Keys_ub_Fmax.begin(), leaf->Keys_ub_Fmax.end(), 0.0);

	leaf->ub_sub = 0.0f;
	leaf->ub_w_in = 0.0f;
	leaf->ub_w_out = 0.0f;

	for (int i = 0; i < leaf->P_s_size; i++) {
		leaf->min_max_dist_P_s[2 * i] = FLT_MAX;
		leaf->min_max_dist_P_s[2 * i + 1] = -FLT_MAX;
	}

	for (User* u : leaf->users) {
		leaf->ub_sub = std::max(leaf->ub_sub, u->ub_sub);
		leaf->ub_w_in = std::max(leaf->ub_w_in, u->ub_w_in);
		leaf->ub_w_out = std::max(leaf->ub_w_out, u->ub_w_out);
		for (int k = 0; k < (int)leaf->Keys_ub_Fsum.size(); ++k) {
			leaf->Keys_ub_Fsum[k] = std::max(leaf->Keys_ub_Fsum[k], u->Keys_Fsum[k]);
			leaf->Keys_ub_Fmax[k] = std::max(leaf->Keys_ub_Fmax[k], u->Keys_Fmax[k]);
		}
		for (int i = 0; i < sd->P_s_size; i++) {
			leaf->min_max_dist_P_s[2 * i] = u->dist_P_s[i] < leaf->min_max_dist_P_s[2 * i] ? u->dist_P_s[i] : leaf->min_max_dist_P_s[2 * i];
			leaf->min_max_dist_P_s[2 * i + 1] = u->dist_P_s[i] > leaf->min_max_dist_P_s[2 * i + 1] ? u->dist_P_s[i] : leaf->min_max_dist_P_s[2 * i + 1];
		}
	}
}
void DynamicApproch::RecomputeInternalMetrics( TreeNode* node)
{
	if (!node) return;

	std::fill(node->Keys_ub_Fsum.begin(), node->Keys_ub_Fsum.end(), 0.0);
	std::fill(node->Keys_ub_Fmax.begin(), node->Keys_ub_Fmax.end(), 0.0);

	node->ub_sub = 0.0f;
	node->ub_w_in = 0.0f;
	node->ub_w_out = 0.0f;
	for (int i = 0; i < node->P_s_size; i++) {
		node->min_max_dist_P_s[2 * i] = FLT_MAX;
		node->min_max_dist_P_s[2 * i + 1] = -FLT_MAX;
	}


	for (TreeNode* ch : node->children) {
		if (!ch) continue;
		node->ub_sub = std::max(node->ub_sub, ch->ub_sub);
		node->ub_w_in = std::max(node->ub_w_in, ch->ub_w_in);
		node->ub_w_out = std::max(node->ub_w_out, ch->ub_w_out);
		for (int k = 0; k < (int)node->Keys_ub_Fsum.size(); ++k) {
			node->Keys_ub_Fsum[k] = std::max(node->Keys_ub_Fsum[k], ch->Keys_ub_Fsum[k]);
			node->Keys_ub_Fmax[k] = std::max(node->Keys_ub_Fmax[k], ch->Keys_ub_Fmax[k]);
		}
		for (int i = 0; i < sd->P_s_size; i++) {
			node->min_max_dist_P_s[2 * i] = ch->min_max_dist_P_s[2 * i] < node->min_max_dist_P_s[2 * i] ? ch->min_max_dist_P_s[2 * i] : node->min_max_dist_P_s[2 * i];
			node->min_max_dist_P_s[2 * i + 1] = ch->min_max_dist_P_s[2 * i + 1] > node->min_max_dist_P_s[2 * i + 1] ? ch->min_max_dist_P_s[2 * i + 1] : node->min_max_dist_P_s[2 * i + 1];
		}
	}
}
void DynamicApproch::UpdateTreeAfterBatch_WithOut_Migration(const std::vector<int>& updated_users)
{
	std::unordered_set<TreeNode*> touched_leaves;

	for (int uid : updated_users) {
		if (uid < 0 || uid >= (int)user_leaf.size()) continue;
		TreeNode* leaf = user_leaf[uid];
		if (leaf) touched_leaves.insert(leaf);
	}

	for (TreeNode* leaf : touched_leaves) {
		RecomputeLeafMetrics(leaf);
	}

	std::unordered_set<TreeNode*> touched_internal;
	for (TreeNode* leaf : touched_leaves) {
		TreeNode* cur = leaf->parent;
		while (cur) {
			touched_internal.insert(cur);
			cur = cur->parent;
		}
	}

	std::vector<TreeNode*> nodes(touched_internal.begin(), touched_internal.end());
	std::sort(nodes.begin(), nodes.end(),
		[](TreeNode* a, TreeNode* b) { return a->depth > b->depth; });

	for (TreeNode* node : nodes) {
		RecomputeInternalMetrics(node);
	}
}

double Normalization_0_1(double x, double mn, double mx)
{
	if (mx <= mn) return 0.0;
	double v = (x - mn) / (mx - mn);
	if (v < 0.0) v = 0.0;
	if (v > 1.0) v = 1.0;
	return v;
}
std::vector<TreeNode*> DynamicApproch::GetAllLeafNodes() const
{
	return leaf_nodes_cached;
}
std::vector<TreeNode*> DynamicApproch::GetNodesAtDepth(int depth) const
{
	auto it = nodes_by_depth.find(depth);
	if (it == nodes_by_depth.end()) return {};
	return it->second;
}
double DynamicApproch::ComputeNodePivotQuality(TreeNode* node, User* piv, float bs_w, float ss_w)
{
	if (!node || !piv || !inx) return 0.0;

	struct NodeRaw7 {
		User* piv = nullptr;
		double bs_sum = 0.0;
		double bs_max = 0.0;
		double ss_sum = 0.0;
		double ss_ub = 0.0;
		double ss_lb = 0.0;
	};

	std::vector<NodeRaw7> scores;
	scores.reserve(P_index_vec.size());

	double mn_rs = DBL_MAX, mn_bs_sum = DBL_MAX, mn_bs_favg = DBL_MAX, mn_ss_sum = DBL_MAX, mn_ss_ub = DBL_MAX, mn_ss_lb = DBL_MAX, mn_bs_max = DBL_MAX;
	double mx_rs = 0.0, mx_bs_sum = 0.0, mx_bs_favg = 0.0, mx_ss_sum = 0.0, mx_ss_ub = 0.0, mx_ss_lb = 0.0, mx_bs_max = 0.0;

	for (User* p : P_index_vec) {
		if (!p) continue;

		NodeRaw7 s;
		s.piv = p;

		auto bs = inx->Get_bs_score_node(node, p);
		s.bs_sum = bs.first;
		s.bs_max = bs.second;

		auto ss = inx->Get_ss_score_node(node, p);
		s.ss_sum = ss.first;
		s.ss_ub = ss.second;
		s.ss_lb = ss.third;

		scores.push_back(s);

		mn_bs_sum = std::min(mn_bs_sum, s.bs_sum);  mx_bs_sum = std::max(mx_bs_sum, s.bs_sum);
		mn_bs_max = std::min(mn_bs_max, s.bs_max);  mx_bs_max = std::max(mx_bs_max, s.bs_max);
		mn_ss_sum = std::min(mn_ss_sum, s.ss_sum);  mx_ss_sum = std::max(mx_ss_sum, s.ss_sum);
		mn_ss_ub = std::min(mn_ss_ub, s.ss_ub);   mx_ss_ub = std::max(mx_ss_ub, s.ss_ub);
		mn_ss_lb = std::min(mn_ss_lb, s.ss_lb);   mx_ss_lb = std::max(mx_ss_lb, s.ss_lb);
	}

	for (const NodeRaw7& s : scores) {
		if (s.piv != piv) continue;

		const double n_bs_sum = Normalization_0_1(s.bs_sum, mn_bs_sum, mx_bs_sum);
		const double n_bs_max = Normalization_0_1(s.bs_max, mn_bs_max, mx_bs_max);
		const double n_ss_sum = Normalization_0_1(s.ss_sum, mn_ss_sum, mx_ss_sum);
		const double n_ss_ub = Normalization_0_1(s.ss_ub, mn_ss_ub, mx_ss_ub);
		const double n_ss_lb = Normalization_0_1(s.ss_lb, mn_ss_lb, mx_ss_lb);

		const double bs_score = 0.5 * (n_bs_sum + n_bs_max);
		const double ss_score = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;
		const double sc = (bs_w * bs_score) + (ss_w * ss_score); // higher is better
		return sc;
	}

	return 0.0;
}

struct UserPivotQualityRow
{
	TreeNode* leaf = nullptr;
	User* piv = nullptr;

	double rs = 0.0;
	double bs_sum = 0.0;
	double bs_favg = 0.0;
	double ss_sum = 0.0;
	double ss_ub = 0.0;
	double ss_lb = 0.0;

	double quality = 0.0;
};
static void NormalizeUserPivotQualities(std::vector<UserPivotQualityRow>& rows, float bs_w, float ss_w, float rs_w)
{
	if (rows.empty()) return;

	double mn_rs = DBL_MAX, mn_bs_sum = DBL_MAX, mn_bs_favg = DBL_MAX;
	double mn_ss_sum = DBL_MAX, mn_ss_ub = DBL_MAX, mn_ss_lb = DBL_MAX;

	double mx_rs = -DBL_MAX, mx_bs_sum = -DBL_MAX, mx_bs_favg = -DBL_MAX;
	double mx_ss_sum = -DBL_MAX, mx_ss_ub = -DBL_MAX, mx_ss_lb = -DBL_MAX;

	for (const auto& r : rows)
	{
		mn_rs = std::min(mn_rs, r.rs);
		mx_rs = std::max(mx_rs, r.rs);

		mn_bs_sum = std::min(mn_bs_sum, r.bs_sum);
		mx_bs_sum = std::max(mx_bs_sum, r.bs_sum);

		mn_bs_favg = std::min(mn_bs_favg, r.bs_favg);
		mx_bs_favg = std::max(mx_bs_favg, r.bs_favg);

		mn_ss_sum = std::min(mn_ss_sum, r.ss_sum);
		mx_ss_sum = std::max(mx_ss_sum, r.ss_sum);

		mn_ss_ub = std::min(mn_ss_ub, r.ss_ub);
		mx_ss_ub = std::max(mx_ss_ub, r.ss_ub);

		mn_ss_lb = std::min(mn_ss_lb, r.ss_lb);
		mx_ss_lb = std::max(mx_ss_lb, r.ss_lb);
	}

	for (auto& r : rows)
	{
		const double n_rs = Normalization_0_1(r.rs, mn_rs, mx_rs);
		const double n_bs_sum = Normalization_0_1(r.bs_sum, mn_bs_sum, mx_bs_sum);
		const double n_bs_favg = Normalization_0_1(r.bs_favg, mn_bs_favg, mx_bs_favg);
		const double n_ss_sum = Normalization_0_1(r.ss_sum, mn_ss_sum, mx_ss_sum);
		const double n_ss_ub = Normalization_0_1(r.ss_ub, mn_ss_ub, mx_ss_ub);
		const double n_ss_lb = Normalization_0_1(r.ss_lb, mn_ss_lb, mx_ss_lb);

		const double bs = 0.5 * (n_bs_sum + n_bs_favg);
		const double ss = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;
		const double rs_sim = 1.0 - n_rs;

		r.quality = (bs_w * bs) + (ss_w * ss) + (rs_w * rs_sim);
	}
}
bool DynamicApproch::MigrateUserToLeafNode(int user_id, float bs_w, float ss_w, float rs_w,	double stability_margin,	const std::vector<TreeNode*>& leaves)
{
	if (!sd || !sd->bn || !inx) return false;
	if (user_id < 0 || user_id >= sd->bn->no_of_users) return false;
	if (user_id >= (int)user_leaf.size()) return false;

	User* u = sd->bn->users_vec[user_id];
	TreeNode* old_leaf = user_leaf[user_id];
	if (!u || !old_leaf) return false;

	if (old_leaf->myUserIndexForPivot < 0 ||
		old_leaf->myUserIndexForPivot >= (int)sd->bn->users_vec.size())
		return false;

	std::vector<UserPivotQualityRow> rows;
	rows.reserve(leaves.size());

	for (TreeNode* leaf : leaves)
	{
		if (!leaf) continue;
		if (leaf->myUserIndexForPivot < 0 ||
			leaf->myUserIndexForPivot >= (int)sd->bn->users_vec.size())
			continue;

		User* piv = sd->bn->users_vec[leaf->myUserIndexForPivot];
		if (!piv) continue;

		UserPivotQualityRow r;
		r.leaf = leaf;
		r.piv = piv;

		r.rs = inx->Get_rs_score(u, piv);

		auto bs = inx->Get_bs_score(u, piv);
		r.bs_sum = bs.first;
		r.bs_favg = bs.second;

		auto ss = inx->Get_ss_score(u, piv);
		r.ss_sum = ss.first;
		r.ss_ub = ss.second;
		r.ss_lb = ss.third;

		rows.push_back(r);
	}

	if (rows.empty()) return false;

	NormalizeUserPivotQualities(rows, bs_w, ss_w, rs_w);

	TreeNode* best_leaf = old_leaf;
	double best_quality = -DBL_MAX;
	double old_quality = -DBL_MAX;

	for (const auto& r : rows)
	{
		if (r.leaf == old_leaf)
			old_quality = r.quality;

		if (r.quality > best_quality)
		{
			best_quality = r.quality;
			best_leaf = r.leaf;
		}
	}

	if (old_quality == -DBL_MAX) return false;
	if (best_leaf == old_leaf) return false;
	if ((best_quality - old_quality) <= stability_margin) return false;


	auto& old_users = old_leaf->users;
	old_users.erase(std::remove(old_users.begin(), old_users.end(), u), old_users.end());


	best_leaf->users.push_back(u);
	user_leaf[user_id] = best_leaf;

	touched_nodes[0].insert(old_leaf);
	touched_nodes[0].insert(best_leaf);

	return true;
}
bool DynamicApproch::MigrateNodeToParentNode( TreeNode* node, float bs_w, float ss_w, double stability_margin)
{
	if (!node || !node->parent) return false;
	if (!sd || !sd->bn) return false;

	TreeNode* old_parent = node->parent;
	const int parent_depth = node->depth - 1;

	std::vector<TreeNode*> parent_candidates = GetNodesAtDepth(parent_depth);
	if (parent_candidates.empty()) return false;

	if (old_parent->myUserIndexForPivot < 0 || old_parent->myUserIndexForPivot >= (int)sd->bn->users_vec.size())
		return false;

	User* old_piv = sd->bn->users_vec[old_parent->myUserIndexForPivot];
	if (!old_piv) return false;

	double old_quality = ComputeNodePivotQuality(node, old_piv, bs_w, ss_w);

	TreeNode* best_parent = old_parent;
	double best_quality = old_quality;


	for (TreeNode* cand_parent : parent_candidates) {
		if (!cand_parent) continue;
		if (cand_parent == old_parent) continue;

		if (cand_parent->depth != parent_depth) continue;
		if (cand_parent->is_leaf) continue;

		if (cand_parent->myUserIndexForPivot < 0 || cand_parent->myUserIndexForPivot >= (int)sd->bn->users_vec.size())
			continue;

		User* piv = sd->bn->users_vec[cand_parent->myUserIndexForPivot];
		if (!piv) continue;

		double quality = ComputeNodePivotQuality(node, piv, bs_w, ss_w);
		if (quality > best_quality) {
			best_quality = quality;
			best_parent = cand_parent;
		}
	}
	if (best_parent == old_parent) return false;
	if ((best_quality - old_quality) <= stability_margin) return false;


	auto& old_children = old_parent->children;
	old_children.erase(std::remove(old_children.begin(), old_children.end(), node), old_children.end());



	auto& new_children = best_parent->children;
	if (std::find(new_children.begin(), new_children.end(), node) == new_children.end()) {
		new_children.push_back(node);
	}
	node->parent = best_parent;
	touched_nodes[old_parent->depth].insert(old_parent);
	touched_nodes[best_parent->depth].insert(best_parent);
	return true;
}
void DynamicApproch::UpdateTreeAfterBatch_With_Migration(float bs_w, float ss_w, float rs_w, double stability_margin, const std::vector<int>& updated_users)
{
	delete inx;
	inx = new Indexing(sd);

	nodes_by_depth.clear();
	leaf_nodes_cached.clear();
	inx = new Indexing(sd);
	for (TreeNode* n : all_nodes) {
		if (!n) continue;
		nodes_by_depth[n->depth].push_back(n);
		if (n->is_leaf)
			leaf_nodes_cached.push_back(n);
	}


	std::vector<TreeNode*> leaves = GetAllLeafNodes();
	P_index_vec.clear();
	for (TreeNode* leaf : leaves) {
		P_index_vec.push_back(sd->bn->users_vec[leaf->myUserIndexForPivot]);
	}

	touched_nodes.clear();
	for (int uid : updated_users) {
		MigrateUserToLeafNode(uid, bs_w, ss_w, rs_w, stability_margin, leaves);
	}

	P_index_vec.clear();
	std::vector<TreeNode*> nodes_dep = GetNodesAtDepth(-1);
	for (TreeNode* nd : nodes_dep) {
		P_index_vec.push_back(sd->bn->users_vec[nd->myUserIndexForPivot]);
	}
	auto it1 = touched_nodes.find(0);
	if (it1 != touched_nodes.end()) {
		std::vector<TreeNode*> cur_nodes(it1->second.begin(), it1->second.end());
		for (TreeNode* n : cur_nodes) {
			RecomputeLeafMetrics(n);
			MigrateNodeToParentNode(n, bs_w, ss_w, stability_margin);
		}
	}

	int test_depth = -1;
	while (true)
	{
		auto it = touched_nodes.find(test_depth);
		if (it == touched_nodes.end() || it->second.empty())
			break;
		std::vector<TreeNode*> cur_nodes(it->second.begin(), it->second.end());

		P_index_vec.clear();
		std::vector<TreeNode*> parent_nodes = GetNodesAtDepth(test_depth - 1);
		for (TreeNode* nd : parent_nodes) {
			P_index_vec.push_back(sd->bn->users_vec[nd->myUserIndexForPivot]);
		}

		for (TreeNode* n : cur_nodes) {
			RecomputeInternalMetrics(n);

			if (n->parent != nullptr) {
				MigrateNodeToParentNode(n, bs_w, ss_w, stability_margin);
			}
		}

		test_depth--;
	}
}



void DynamicApproch::RebuildUserWindowedCheckins(User* u, int current_day, int window_days)
{
	u->checkin_locations.clear();

	const int expire_before = current_day - window_days;

	for (const auto& kv : u->timed_checkins) {
		int poi_id = kv.first;
		const std::vector<int>& days = kv.second;

		for (int day : days) {
			if (expire_before <= day && day <= current_day) {
				u->checkin_locations[poi_id]++;
			}
		}
	}
}
void DynamicApproch::RebuildUserKeyStats(User* u)
{
	if (!sd) return;

	const int K = sd->toatl_n_of_used_keys;

	u->Keys_Fsum.assign(K, 0.0);
	u->Keys_Fmax.assign(K, 0.0);
	u->Keys_visited_count.assign(K, 0);

	for (const auto& pf : u->checkin_locations) {
		int poi_id = pf.first;
		int f = pf.second;
		if (f <= 0) continue;

		POI* p = sd->poi_vec[poi_id];
		if (!p) continue;

		for (int i = 0; i < p->num_of_keys; ++i) {
			int key = p->keys[i];
			if (key < 0 || key >= K) continue;

			u->Keys_Fsum[key] += f;
			u->Keys_visited_count[key] += 1;
			u->Keys_Fmax[key] = std::max(u->Keys_Fmax[key], (double)f);
		}
	}
}
void DynamicApproch::RecomputeUserUb(User* u)
{
	float ub_avg = 0.0f;
	int ub_sum = 0;

	for (const auto& pf : u->checkin_locations) {
		ub_avg = std::max(ub_avg, (float)pf.second);
		ub_sum += pf.second;
	}

	u->ub_f_avg = ub_avg;
	u->ub_f_sum = ub_sum;
}

void DynamicApproch::GenerateOfflineData1( int pivot_iterations, std::string filename, int distribution_type) {

		double spent_time = 0.0;

	OfflineCalculations* offlineCal = new OfflineCalculations(sd);
	offlineCal->Run(pivot_iterations, filename, distribution_type);
	delete offlineCal;
}
std::vector<CommStru> DynamicApproch::Load_Communities(std::string filename) {
	std::vector<CommStru> comm_vec;


	std::string line;
	std::ifstream fileComm("dataset/comm/" + filename + "_communities.txt");

	if (!fileComm.is_open())
	{
		std::cout << "Could not open file " << filename << " !\n";
		return comm_vec;
	}
	std::string users_str, poi_str;
	
	while (getline(fileComm, line))
	{
		if (line.length() < 1) continue;
		CommStru C;
		C.q_index = std::stoi(line.substr(0, line.find(":")));
		line.erase(0, line.find(":") + 1);
		users_str = line.substr(0, line.find("&"));
		line.erase(0, line.find("&") + 1);
		poi_str = line.substr(0, line.find(" "));
		line.erase(0, line.find(" ") + 1);


		while (users_str.length() > 1) {
			C.user_ids.insert(std::stoi(users_str.substr(0, users_str.find(","))));
			users_str.erase(0, users_str.find(",") + 1);
		}

		while (poi_str.length() > 1) {
			C.poi_ids.insert(std::stoi(poi_str.substr(0, poi_str.find(","))));
			poi_str.erase(0, poi_str.find(",") + 1);
		}
		comm_vec.push_back(C);
	}
	fileComm.close();

	return comm_vec;
}
void DynamicApproch::Add_random_timestamps_for_first_time( std::string filename) {
	std::unordered_map<int, std::vector<int>> timed_checkins;
	std::vector<int> timestamps;
	for (int u_id = 0; u_id < sd->bn->no_of_users; u_id++) {
		int checkin_count = sd->bn->users_vec[u_id]->checkin_locations.size();
		for (auto p_f : sd->bn->users_vec[u_id]->checkin_locations) {
			for (int i = 0; i < p_f.second; i++) {
				timestamps.push_back(random1::uniform_int_pos(1, 50));
			}
			sort(timestamps.begin(), timestamps.end());
			sd->bn->users_vec[u_id]->timed_checkins[p_f.first] = timestamps;
			timestamps.clear();
		}
	}
	sd->Print_SN_timestamp(true, filename);
}
std::vector<UserVisit> DynamicApproch::Picked_users_get_visits( std::vector<int> picked_users, int max_visits, int min_day, int max_day) {
	std::vector<UserVisit> user_visits;

	for (int i = 0; i < picked_users.size(); i++) {
		UserVisit uv;
		uv.user_id = picked_users[i];
		uv.poi_id = random1::uniform_int_pos(0, sd->no_of_POI - 1); 
		int random_visit_count = random1::uniform_int_pos(1, max_visits); 
		uv.days.reserve(random_visit_count);
		for (int j = 0; j < random_visit_count; j++) {
			uv.days.push_back(random1::uniform_int_pos(min_day, max_day));
		}
		user_visits.push_back(uv);

	}
	return user_visits;
}
void DynamicApproch::Picked_users_Add_visit( const std::vector<UserVisit>& user_visits) {
	for (int i = 0; i < user_visits.size(); i++) {
		User* user = sd->bn->users_vec[user_visits[i].user_id];
		user->checkin_locations[user_visits[i].poi_id] += user_visits[i].days.size(); 		
	}
}
std::vector<int> DynamicApproch::Picked_users_Filter_visit( std::vector<int> picked_users, int max_day, int window_days) {

	std::vector<int> updated_users;
	for (int u_id : picked_users) {
		User* u = sd->bn->users_vec[u_id];
		std::map<int, int> old_checkins = u->checkin_locations;
		RebuildUserWindowedCheckins(u, max_day, window_days);
		if (u->checkin_locations != old_checkins)
			updated_users.push_back(u->id);
	}
	return updated_users;
}
void DynamicApproch::Reset_Remaining_By_Community( CommStru& comm) {

	sd->no_remining_users = 0;
	sd->no_remining_POI = 0;
	for (int i = 0; i < sd->bn->no_of_users; i++)
		sd->bn->users_vec[i]->isPruned = true;
	for (int i = 0; i < sd->no_of_POI; i++)
		sd->poi_vec[i]->isPruned = true;
	for (int u_id : comm.user_ids) {
		sd->bn->users_vec[u_id]->isPruned = false;
		sd->remining_users[sd->no_remining_users] = u_id;
		sd->no_remining_users++;
	}
	for (int p_id : comm.poi_ids) {
		sd->poi_vec[p_id]->isPruned = false;
		sd->remining_POI[sd->no_remining_POI] = p_id;
		sd->no_remining_POI++;
	}
}
void DynamicApproch::Add_affected_Users_With_POI_to_Remaining_all(std::vector<int> affected_users, CommStru& comm) {

	for (int u_id : affected_users) {
		if (comm.user_ids.count(u_id) == 0) {
			if (sd->bn->users_vec[u_id]->isPruned) {
				sd->bn->users_vec[u_id]->isPruned = false;
				sd->remining_users[sd->no_remining_users] = u_id;
				sd->no_remining_users++;
			}
			comm.user_ids.insert(u_id);
		}
		for (auto p_f : sd->bn->users_vec[u_id]->checkin_locations) {
			if (comm.poi_ids.count(p_f.first) == 0) {
				int p_id = p_f.first;
				if (sd->poi_vec[p_id]->isPruned) {
					sd->poi_vec[p_id]->isPruned = false;
					sd->remining_POI[sd->no_remining_POI] = p_id;
					sd->no_remining_POI++;
				}
			}
			comm.poi_ids.insert(p_f.first);
		}
	}
}
void DynamicApproch::Add_affected_Users_With_POI_to_Remaining_one( int u_id, CommStru& comm) {

	if (sd->bn->users_vec[u_id]->isPruned) {
		sd->bn->users_vec[u_id]->isPruned = false;
		sd->remining_users[sd->no_remining_users] = u_id;
		sd->no_remining_users++;
	}
	comm.user_ids.insert(u_id);
	for (auto p_f : sd->bn->users_vec[u_id]->checkin_locations) {
		if (comm.poi_ids.count(p_f.first)) continue;
		comm.poi_ids.insert(p_f.first);
		int p_id = p_f.first;
		if (sd->poi_vec[p_id]->isPruned) {
			sd->poi_vec[p_id]->isPruned = false;
			sd->remining_POI[sd->no_remining_POI] = p_id;
			sd->no_remining_POI++;
		}
	}

}
bool DynamicApproch::Validate_Communities_After_Repair(std::string filename, Ref2* ref, int q_index,
	std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value) {

	int omega_int = omega_value * sd->larget_f_sum;
	float pi_fl = pi_value * sd->larget_f_avg;
	if (pi_fl < 1)pi_fl = 1;
	ref->Run_Ref_dynamic(Q, k_value, d_value, omega_int, pi_fl, theta_value, sigma_value, q_index);
	return sd->no_remining_users > 2;
}
std::unordered_set<int> DynamicApproch::Select_random_comm_index_to_test(std::string filename, int num_to_select) {

	std::vector<CommStru> comm_vec = Load_Communities(filename);
	std::unordered_set<int> random_comm_index_to_test_set;

	if (comm_vec.size() > num_to_select)
		for (int i = 0; i < num_to_select; i++) {
			int index = random1::uniform_int_pos(0, comm_vec.size() - 1);
			while (random_comm_index_to_test_set.count(index) != 0)
				index = random1::uniform_int_pos(0, comm_vec.size() - 1);
			random_comm_index_to_test_set.insert(index);
		}
	else {
		for (int i = 0; i < comm_vec.size(); i++)
			random_comm_index_to_test_set.insert(i);
	}
	return random_comm_index_to_test_set;
}
void DynamicApproch::Build_tree_and_save_withComm( std::string filename, bool loadTree,
	std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value) {

	int pivot_iterations = 5;
	GenerateOfflineData1(pivot_iterations, filename, 1);
	int tree_piv_size[] = { 100 };
	Tree* tree = NULL;
	if (loadTree)
		tree = TreeTextIO::Load(sd,"dataset/tree/"+ filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	else {
		float bs_w = 1, ss_w = 2, rs_w = 1;
		int tree_piv_size_single = sd->bn->no_of_users / 100, level_dev = 2;
		tree = new Tree(sd);
		tree->BuildTree(pivot_iterations, level_dev, tree_piv_size_single, bs_w, ss_w, rs_w);
		TreeTextIO::Save(tree, filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	}

	Ref2* ref = new Ref2(sd);
	
	FindCommuntiesPrintToFile_Defults(tree, filename, ref, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value);

	delete tree;
	delete ref;
}
std::string DynamicApproch::FindCommuntiesPrintToFile_Defults( Tree* tree, std::string filename, Ref2* ref,
	std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value) {



	int omega_int = omega_value * sd->larget_f_sum;
	float pi_fl = pi_value * sd->larget_f_avg;
	if (pi_fl < 1)pi_fl = 1;

	int q_tested = 0;
	int text_all = 0;
	int q_index;
	std::ostringstream comm;

	q_index = -1;
	while (q_index < sd->bn->no_of_users - 1) {
		q_index++;

		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		sd->Reset_POI_pruned_flag();
		sd->bn->Reset_Users_pruned_flag();
		for (auto u_id : sd->bn->users_vec) {
			sd->remining_users[sd->no_remining_users] = u_id->id;
			sd->no_remining_users++;
		}
		int remain_users2 = 0;
		int remian_poi = 0;

		tree->KCS_BSSN_Query_Answer(Q, k_value, d_value, omega_int, pi_fl, theta_value, sigma_value, q_index, { true,true, true, true, true, true, true });
		ref->Run_Ref2(Q, k_value, d_value, omega_int, pi_fl, theta_value, sigma_value, q_index, { true,true, true, true, true, true, true });

		remain_users2 = sd->no_remining_users;
		remian_poi = sd->no_remining_POI;
		if (sd->no_remining_users < 3) continue;



		q_tested++;
		comm << q_index << ":";
		for (int i = 0; i < sd->no_remining_users; i++) {
			comm << sd->remining_users[i] << ",";
		}
		comm << "&";
		for (int i = 0; i < sd->no_remining_POI; i++) {
			comm << sd->remining_POI[i] << ",";
		}
		comm << "\n";

	}
	
	sd->Print_To_File(filename + "_communities.txt", comm.str());
	return "";

}
std::string DynamicApproch::Run_Dynamic_Temporal_Insert_visits(SampleData* _sd, int user_no_updates, std::vector<int> updated_users, int max_visits, int min_day, int max_day, std::string filename,
	std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value,
	float bs_w, float ss_w, float rs_w, double stability_margin, std::unordered_set<int> comm_indeies) {
	sd = _sd;

	Ref2* ref = new Ref2(sd);
	double repair_dataset_time = 0.0;
	const auto start = std::chrono::high_resolution_clock::now();


	std::vector<UserVisit> user_visits = Picked_users_get_visits( updated_users, max_visits, min_day, max_day);


	
	Picked_users_Add_visit(user_visits);
	const auto end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms = end - start;
	repair_dataset_time = duration_ms.count();


	int pivot_iterations = 5;
	GenerateOfflineData1( pivot_iterations, filename, 1);
	


	double offline_time = 0.0;
	const auto start1 = std::chrono::high_resolution_clock::now();
	for (int user_id : updated_users) {
		User* u = sd->bn->users_vec[user_id];
		RebuildUserKeyStats(u);
		RecomputeUserUb(u);
	}
	const auto end1 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms1 = end1 - start1;
	offline_time = duration_ms1.count();




	double Comm_load_time = 0.0;
	const auto start2 = std::chrono::high_resolution_clock::now();
	std::vector<CommStru> comm_vec = Load_Communities(filename);
	const auto end2 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms2 = end2 - start2;
	Comm_load_time = duration_ms2.count();

	std::vector<CommStru> comm_vec_copy = comm_vec;

	
	double Comm_repair_time_total = 0.0;
	const auto start3 = std::chrono::high_resolution_clock::now();
	for (auto comm_index : comm_indeies) {
		if (comm_index < 0 || comm_index >= (int)comm_vec.size()) continue;
		Reset_Remaining_By_Community( comm_vec[comm_index]);
		for (int u_id : updated_users) {
			bool user_in_comm = comm_vec[comm_index].user_ids.count(u_id) != 0;
			if (user_in_comm) continue; 
			Add_affected_Users_With_POI_to_Remaining_one( u_id, comm_vec[comm_index]);
			Validate_Communities_After_Repair( filename, ref, comm_vec[comm_index].q_index, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value);
		}
	}
	const auto end3 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms3 = end3 - start3;
	Comm_repair_time_total = duration_ms3.count();
	double Comm_repair_time_avg = 0.0;
	if (!comm_indeies.empty())
		Comm_repair_time_avg = Comm_repair_time_total / comm_indeies.size();


	
	double Comm_repair_time_batch_total = 0.0;
	const auto start31 = std::chrono::high_resolution_clock::now();
	for (auto comm_index : comm_indeies) {
		if (comm_index < 0 || comm_index >= (int)comm_vec.size()) continue;
		Reset_Remaining_By_Community( comm_vec_copy[comm_index]);
		Add_affected_Users_With_POI_to_Remaining_all( updated_users, comm_vec_copy[comm_index]);
		Validate_Communities_After_Repair(  filename, ref, comm_vec_copy[comm_index].q_index, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value);
	}
	const auto end31 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms31 = end31 - start31;
	Comm_repair_time_batch_total = duration_ms31.count();
	double Comm_repair_time_batch_avg = 0.0;
	if (!comm_indeies.empty())
		Comm_repair_time_batch_avg = Comm_repair_time_batch_total / comm_indeies.size();



	Tree* tree = NULL;
	tree = TreeTextIO::Load(sd, filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	all_nodes = tree->all_nodes;
	user_leaf = tree->user_leaf;


	
	double tree_node_repair_time = 0.0;
	const auto start4 = std::chrono::high_resolution_clock::now();
	UpdateTreeAfterBatch_WithOut_Migration(updated_users);
	const auto end4 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms4 = end4 - start4;
	tree_node_repair_time = duration_ms4.count();

	delete tree;
	tree = NULL;
	tree = TreeTextIO::Load(sd, filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	all_nodes = tree->all_nodes;
	user_leaf = tree->user_leaf;


	double tree_node_repair_with_migrate_time = 0.0;
	const auto start5 = std::chrono::high_resolution_clock::now();
	UpdateTreeAfterBatch_With_Migration(bs_w, ss_w, rs_w, stability_margin, updated_users);
	const auto end5 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms5 = end5 - start5;
	tree_node_repair_with_migrate_time = duration_ms5.count();


	delete tree;
	delete ref;
	
	
	return  (std::to_string(user_no_updates) + "," + std::to_string(updated_users.size()) + "," + std::to_string(repair_dataset_time) + "," + std::to_string(offline_time) + "," + std::to_string(Comm_load_time) + "," + std::to_string(comm_indeies.size()) + "," + std::to_string(Comm_repair_time_total) + "," + std::to_string(Comm_repair_time_avg) + "," + std::to_string(Comm_repair_time_batch_total) + "," + std::to_string(Comm_repair_time_batch_avg) + "," + std::to_string(tree_node_repair_time) + "," + std::to_string(tree_node_repair_with_migrate_time));

}

std::string DynamicApproch::Run_Dynamic_Temporal_Delete_visits_by_user_number(SampleData* _sd, int user_no_updates, std::vector<int> picked_users, int max_day, int window_days, std::string filename,
	std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value,
	float theta_value, float sigma_value, float bs_w, float ss_w, float rs_w, double stability_margin, std::unordered_set<int> comm_indeies) {
	sd = _sd;

	Ref2* ref = new Ref2(sd);
	double repair_dataset_time = 0.0;
	const auto start = std::chrono::high_resolution_clock::now();

	
	std::vector<int> updated_users;
	updated_users = Picked_users_Filter_visit( picked_users, max_day, window_days);
	const auto end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms = end - start;
	repair_dataset_time = duration_ms.count();


	int pivot_iterations = 5;
	GenerateOfflineData1( pivot_iterations, filename, 1);
	


	double offline_time = 0.0;
	const auto start1 = std::chrono::high_resolution_clock::now();
	for (int user_id : updated_users) {
		User* u = sd->bn->users_vec[user_id];
		RebuildUserKeyStats(u);
		RecomputeUserUb(u);
	}
	const auto end1 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms1 = end1 - start1;
	offline_time = duration_ms1.count();


	
	double Comm_load_time = 0.0;
	const auto start2 = std::chrono::high_resolution_clock::now();
	std::vector<CommStru> comm_vec = Load_Communities(filename);
	const auto end2 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms2 = end2 - start2;
	Comm_load_time = duration_ms2.count();

	std::vector<CommStru> comm_vec_copy = comm_vec;




	
	double Comm_repair_time_total = 0.0;
	const auto start3 = std::chrono::high_resolution_clock::now();
	for (auto comm_index : comm_indeies) {
		if (comm_index < 0 || comm_index >= (int)comm_vec.size()) continue;
		Reset_Remaining_By_Community( comm_vec[comm_index]);
		for (int u_id : updated_users) {
			if (comm_vec[comm_index].user_ids.count(u_id)==0) continue; 
			Add_affected_Users_With_POI_to_Remaining_one( u_id, comm_vec[comm_index]);
			Validate_Communities_After_Repair(  filename, ref, comm_vec[comm_index].q_index, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value);
		}


	}
	const auto end3 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms3 = end3 - start3;
	Comm_repair_time_total = duration_ms3.count();
	double Comm_repair_time_avg = 0.0;
	if (!comm_indeies.empty())
		Comm_repair_time_avg = Comm_repair_time_total / comm_indeies.size();


	
	double Comm_repair_time_batch_total = 0.0; const auto start31 = std::chrono::high_resolution_clock::now();
	for (auto comm_index : comm_indeies) {
		if (comm_index < 0 || comm_index >= (int)comm_vec.size()) continue;
		Reset_Remaining_By_Community( comm_vec_copy[comm_index]);
		Add_affected_Users_With_POI_to_Remaining_all( updated_users, comm_vec_copy[comm_index]);
		Validate_Communities_After_Repair(  filename, ref, comm_vec_copy[comm_index].q_index, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value);

	}
	const auto end31 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms31 = end31 - start31;
	Comm_repair_time_batch_total = duration_ms31.count();
	double Comm_repair_time_batch_avg = 0.0;
	if (!comm_indeies.empty())
		Comm_repair_time_batch_avg = Comm_repair_time_batch_total / comm_indeies.size();

	Tree* tree = NULL;
	tree = TreeTextIO::Load(sd, filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	all_nodes = tree->all_nodes;
	user_leaf = tree->user_leaf;



	double tree_node_repair_time = 0.0;
	const auto start4 = std::chrono::high_resolution_clock::now();
	UpdateTreeAfterBatch_WithOut_Migration(updated_users);
	const auto end4 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms4 = end4 - start4;
	tree_node_repair_time = duration_ms4.count();

	delete tree;
	tree = NULL;
	tree = TreeTextIO::Load(sd, filename + "_" + std::to_string(sd->bn->no_of_users) + "_tree_cache.txt");
	all_nodes = tree->all_nodes;
	user_leaf = tree->user_leaf;


	
	double tree_node_repair_with_migrate_time = 0.0;
	const auto start5 = std::chrono::high_resolution_clock::now();
	UpdateTreeAfterBatch_With_Migration(bs_w, ss_w, rs_w, stability_margin, updated_users);
	const auto end5 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms5 = end5 - start5;
	tree_node_repair_with_migrate_time = duration_ms5.count();


	delete tree;
	delete ref;
	return (std::to_string(user_no_updates) + "," + std::to_string(updated_users.size()) + "," + std::to_string(repair_dataset_time) + "," + std::to_string(offline_time) + "," + std::to_string(Comm_load_time) + "," + std::to_string(comm_indeies.size()) + "," + std::to_string(Comm_repair_time_total) + "," + std::to_string(Comm_repair_time_avg) + "," + std::to_string(Comm_repair_time_batch_total) + "," + std::to_string(Comm_repair_time_batch_avg) + "," + std::to_string(tree_node_repair_time) + "," + std::to_string(tree_node_repair_with_migrate_time));

}











