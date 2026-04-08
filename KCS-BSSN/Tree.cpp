#include "Tree.h"
#include "rand.h"
#include <queue>
#include "Pruning.h"
#include "MinHeap.h"
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <limits>

Tree::Tree(SampleData* _sd) {
	sd = _sd;
	num_of_nodes = 0;
	root = NULL;

	if (sd && sd->bn) {
		user_leaf.resize(sd->bn->no_of_users, nullptr);
	}

}
Tree::~Tree()
{
	for (TreeNode* n : all_nodes) {
		delete n;
	}
	all_nodes.clear();
	Nodes.clear();
	root = NULL;
}
void Tree::BulidLeafNodes() {
	int i = 1;
	for (std::pair < User*, std::vector<User*>> group : inx->Subgraphs) {
		TreeNode* tnode = Create_Leaf_Node(inx, group.first->id, group.second);
		tnode->id = i++;
		Nodes.push_back(tnode);
		num_of_nodes = num_of_nodes + 1;
	}
}
void Tree::BuildTree(int iteration, int  tree_level_division, int piv_size, float bs_w, float ss_w, float rs_w)
{
	inx = new Indexing(sd);
	inx->Indexing_Pivot_Selection(iteration, piv_size, bs_w, ss_w, rs_w);
	for (int i = 0; i < inx->P_index_vec.size(); i++)
		P_index_vec_for_nodes.push_back(inx->P_index_vec[i]);
	BulidLeafNodes();
	ButtomUp(tree_level_division, iteration, bs_w, ss_w, rs_w);

	if (Nodes.empty()) {
		root = NULL;
	}
	else if (Nodes.size() == 1) {
		root = Nodes[0];
		root->isroot = true;
	}
	else {
		TreeNode* super_root = Create_Non_Leaf_Node(Nodes[0]->depth - 1, Nodes, -1);
		Nodes.clear();
		Nodes.push_back(super_root);
		++num_of_nodes;
		root = super_root;
		root->isroot = true;
	}
	delete inx;
	inx = nullptr;
}

TreeNode* Tree::Create_Leaf_Node(Indexing* inx, int pivot_user_index, const std::vector<User*>& group)
{
	TreeNode* tnode_new = new TreeNode(0, sd->P_s_size, sd->P_r_size, this, true, pivot_user_index, sd->toatl_n_of_used_keys);
	all_nodes.push_back(tnode_new);
	tnode_new->id = all_nodes.size();
	for (User* u : group) {
		tnode_new->users.push_back(u);
		user_leaf[u->id] = tnode_new;

		tnode_new->ub_sub = tnode_new->ub_sub < u->ub_sub ? u->ub_sub : tnode_new->ub_sub;
		tnode_new->ub_w_in = u->ub_w_in > tnode_new->ub_w_in ? u->ub_w_in : tnode_new->ub_w_in;
		tnode_new->ub_w_out = u->ub_w_out > tnode_new->ub_w_out ? u->ub_w_out : tnode_new->ub_w_out;
		for (int i = 0; i < sd->toatl_n_of_used_keys; i++) {
			tnode_new->Keys_ub_Fsum[i] = fmax(tnode_new->Keys_ub_Fsum[i], u->Keys_Fsum[i]);
			tnode_new->Keys_ub_Fmax[i] = fmax(tnode_new->Keys_ub_Fmax[i], u->Keys_Fmax[i]);
		}
		for (int i = 0; i < sd->P_s_size; i++) {
			tnode_new->min_max_dist_P_s[2 * i] = u->dist_P_s[i] < tnode_new->min_max_dist_P_s[2 * i] ? u->dist_P_s[i] : tnode_new->min_max_dist_P_s[2 * i];
			tnode_new->min_max_dist_P_s[2 * i + 1] = u->dist_P_s[i] > tnode_new->min_max_dist_P_s[2 * i + 1] ? u->dist_P_s[i] : tnode_new->min_max_dist_P_s[2 * i + 1];
		}

	}
	return tnode_new;
}
TreeNode* Tree::Create_Non_Leaf_Node(int level, std::vector<TreeNode*> Children, int pivot_user_index)
{
	TreeNode* tnode_new = new TreeNode(level, sd->P_s_size, sd->P_r_size, this, false, pivot_user_index, sd->toatl_n_of_used_keys);
	all_nodes.push_back(tnode_new);
	tnode_new->id = all_nodes.size();
	tnode_new->children = Children;
	for (TreeNode* child : Children) {
		child->parent = tnode_new;
		tnode_new->ub_sub = tnode_new->ub_sub < child->ub_sub ? child->ub_sub : tnode_new->ub_sub;
		tnode_new->ub_w_in = tnode_new->ub_w_in < child->ub_w_in ? child->ub_w_in : tnode_new->ub_w_in;
		tnode_new->ub_w_out = tnode_new->ub_w_out < child->ub_w_out ? child->ub_w_out : tnode_new->ub_w_out;
		for (int i = 0; i < sd->toatl_n_of_used_keys; i++) {
			tnode_new->Keys_ub_Fsum[i] = fmax(tnode_new->Keys_ub_Fsum[i], child->Keys_ub_Fsum[i]);
			tnode_new->Keys_ub_Fmax[i] = fmax(tnode_new->Keys_ub_Fmax[i], child->Keys_ub_Fmax[i]);
		}
		for (int i = 0; i < sd->P_s_size; i++) {
			tnode_new->min_max_dist_P_s[2 * i] = child->min_max_dist_P_s[2 * i] < tnode_new->min_max_dist_P_s[2 * i] ? child->min_max_dist_P_s[2 * i] : tnode_new->min_max_dist_P_s[2 * i];
			tnode_new->min_max_dist_P_s[2 * i + 1] = child->min_max_dist_P_s[2 * i + 1] > tnode_new->min_max_dist_P_s[2 * i + 1] ? child->min_max_dist_P_s[2 * i + 1] : tnode_new->min_max_dist_P_s[2 * i + 1];
		}
	}
	return tnode_new;
}

void Tree::ButtomUp(int tree_level_division, int piv_selction_iter, float bs_w, float ss_w, float rs_w)
{
	const size_t prev_size = Nodes.size();
	if (prev_size <= 1) return;

	int piv_size = (int)prev_size / tree_level_division;
	if (piv_size < 1) piv_size = 1;

	std::vector<User*> new_p_index_set =
		inx->Indexing_Pivot_Selection_node(piv_selction_iter, Nodes, P_index_vec_for_nodes, piv_size, bs_w, ss_w);

	P_index_vec_for_nodes = new_p_index_set;

	struct Scores { double bs_sum, bs_max, ss_sum, ss_ub, ss_lb; };
	std::vector<Scores> scores;
	scores.resize(P_index_vec_for_nodes.size());

	std::map<int, std::vector<TreeNode*>> Piv_Nodes_map;
	for (TreeNode* N : Nodes) {

		Scores maxes{ 0,0,0,0,0 };
		Scores mins{ DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX };


		for (size_t t = 0; t < P_index_vec_for_nodes.size(); ++t) {
			User* piv = P_index_vec_for_nodes[t];

			auto fsum_fmax = inx->Get_bs_score_node(N, piv);
			auto ss_score = inx->Get_ss_score_node(N, piv);

			Scores s{ fsum_fmax.first, fsum_fmax.second, ss_score.first, ss_score.second, ss_score.third };
			scores[t] = s;

			maxes.bs_sum = std::max(maxes.bs_sum, s.bs_sum);
			maxes.bs_max = std::max(maxes.bs_max, s.bs_max);
			maxes.ss_sum = std::max(maxes.ss_sum, s.ss_sum);
			maxes.ss_ub = std::max(maxes.ss_ub, s.ss_ub);
			maxes.ss_lb = std::max(maxes.ss_lb, s.ss_lb);

			mins.bs_sum = std::min(mins.bs_sum, s.bs_sum);
			mins.bs_max = std::min(mins.bs_max, s.bs_max);
			mins.ss_sum = std::min(mins.ss_sum, s.ss_sum);
			mins.ss_ub = std::min(mins.ss_ub, s.ss_ub);
			mins.ss_lb = std::min(mins.ss_lb, s.ss_lb);
		}


		double best_quality = -1.0;
		int best_pivot_id = P_index_vec_for_nodes[0]->id;

		auto norm = [](double x, double mn, double mx) -> double {
			return (mx > mn) ? ((x - mn) / (mx - mn)) : 0.0;
			};

		for (size_t t = 0; t < P_index_vec_for_nodes.size(); ++t) {
			const Scores& s = scores[t];

			const double n_bs_sum = norm(s.bs_sum, mins.bs_sum, maxes.bs_sum);
			const double n_bs_max = norm(s.bs_max, mins.bs_max, maxes.bs_max);
			const double n_ss_sum = norm(s.ss_sum, mins.ss_sum, maxes.ss_sum);
			const double n_ss_ub = norm(s.ss_ub, mins.ss_ub, maxes.ss_ub);
			const double n_ss_lb = norm(s.ss_lb, mins.ss_lb, maxes.ss_lb);

			const double bs_score = (n_bs_sum + n_bs_max) / 2.0;
			const double ss_score = (n_ss_sum + n_ss_ub + (1.0 - n_ss_lb)) / 3.0;
			const double quality = ((bs_w * bs_score) + (ss_w * ss_score)) / 2.0;

			if (quality > best_quality) {
				best_quality = quality;
				best_pivot_id = P_index_vec_for_nodes[t]->id;
			}
		}

		Piv_Nodes_map[best_pivot_id].push_back(N);
	}

	Nodes.clear();
	Nodes.reserve(Piv_Nodes_map.size());

	for (auto& kv : Piv_Nodes_map) {
		auto& group = kv.second;
		TreeNode* non_leaf_node = Create_Non_Leaf_Node(group[0]->depth - 1, group, kv.first);
		Nodes.push_back(non_leaf_node);
		++num_of_nodes;
	}

	if (Nodes.size() >= prev_size) return; 
	ButtomUp(tree_level_division, piv_selction_iter, bs_w, ss_w, rs_w);
}

bool Tree::KCS_BSSN_Query_Answer(std::vector<int> Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index, Sextuple<bool, bool, bool, bool, bool, bool, bool> m)
{
	for (auto* u : sd->bn->users_vec) if (u) u->isPruned = true;
	for (auto* p : sd->poi_vec)       if (p) p->isPruned = false;
	sd->no_remining_users = 0;
	sd->no_remining_POI = 0;

	Pruning pruning(sd);

	User* q = sd->bn->users_vec[q_index];
	bool keyword_pruning = m.keyword && pruning.Keyword_based_pruning(q, Q);
	bool omega_pruning = m.omega && pruning.Omega_based_pruning(q, omega);
	bool pi_pruning = m.pi && pruning.Pi_based_pruning(q, pi);
	bool structural_cohesiveness_pruning = m.structrual && pruning.Structural_cohesiveness_pruning_sub(q, k);

	if (keyword_pruning || omega_pruning || pi_pruning || structural_cohesiveness_pruning) {
		sd->no_remining_POI = 0;
		sd->no_remining_users = 0;
		return true;
	}
	std::queue<TreeNode*> queueNode;
	std::queue<User*> queueUser;

	if (!root) {
		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		return true;
	}
	queueNode.push(root);

	int i = 0;
	while (!queueNode.empty()) {
		TreeNode* node = queueNode.front();
		queueNode.pop();

		bool keyword_pruning_node = m.keyword && pruning.Keyword_based_pruning_index_node(node, Q, d);
		bool pi_pruning_node = m.pi && pruning.Pi_based_pruning_index_node(node, Q, pi, d);
		bool omega_pruning_node = m.omega && pruning.Omega_based_pruning_index_node(node, Q, omega, d);
		bool structural_cohesiveness_pruning_node = m.structrual && pruning.Structural_cohesiveness_pruning_sub_index_node(node, k);
		bool influence_pruning_node = m.influnce && pruning.Influnce_based_pruning_index_node(node, q, theta);
		bool social_distance_pruning_node = m.social && pruning.Social_distance_based_pruning_index_node(node, q, d);
		if (keyword_pruning_node || omega_pruning_node || pi_pruning_node || influence_pruning_node || structural_cohesiveness_pruning_node || social_distance_pruning_node)
			continue;


		if (!node->is_leaf)
			for (TreeNode* child : node->children)
				queueNode.push(child);
		else
			for (User* child_user : node->users)
			{
				if (child_user->id != q->id) {
					bool keyword_pruning = m.keyword && pruning.Keyword_based_pruning(child_user, Q);
					bool pi_pruning = m.pi && pruning.Pi_based_pruning(child_user, pi);
					bool omega_pruning = m.omega && pruning.Omega_based_pruning(child_user, omega);
					bool spatial_distance_pruning = m.spacial && pruning.Spatial_distance_based_pruning(child_user, q, Q, sigma);
					bool structural_cohesiveness_pruning = m.structrual && pruning.Structural_cohesiveness_pruning_sub(child_user, k);
					bool influence_pruning = m.influnce && pruning.Influnce_based_pruning(q, child_user, theta);
					bool social_distance_pruning = m.social && pruning.Social_distance_based_pruning(child_user, q, d);
					if (keyword_pruning || omega_pruning || pi_pruning || influence_pruning || structural_cohesiveness_pruning || social_distance_pruning || spatial_distance_pruning)
					{

						continue;
					}
				}
				sd->remining_users[i] = child_user->id;
				sd->bn->users_vec[child_user->id]->isPruned = false;
				i++;

			}
	}
	sd->no_remining_users = i;
	return false;
}




