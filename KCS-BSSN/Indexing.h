#pragma once
#include <vector>
#include "User.h"
#include "SampleData.h"
#include "TreeNode.h"
#include "Triple.h"
#include "Sextuple.h"
class TreeNode;
class Indexing
{
public:

	SampleData* sd;
	std::vector<User*> P_index_vec;
	std::map < User*, std::vector<User*>> Subgraphs;// pivot_index, subgraph	
	Indexing(SampleData* _sd);
	~Indexing();

	std::pair<double, double> Get_bs_score(User* u, User* v);
	double Get_rs_score(User* u, User* v);
	Triple<double, double, double> Get_ss_score(User* u, User* v);
	std::pair<double, double> Get_bs_score_node(TreeNode* N, User* v);
	Triple<double, double, double> Get_ss_score_node(TreeNode* N, User* v);

	void Indexing_Pivot_Selection(int iteration, int piv_size, float bs_w, float ss_w, float rs_w);
	void Partition_SN(float bs_w, float ss_w, float rs_w);
	double calculate_p_index_cost(float bs_w, float ss_w, float rs_w);
	std::vector<User*> Indexing_Pivot_Selection_node(int iteration, const std::vector<TreeNode*>& Nodes, const std::vector<User*>& P_index_vec, int piv_size, float bs_w, float ss_w);
	float calculate_p_index_cost_node(const std::vector<TreeNode*>& Nodes, const std::vector<User*>& temp_P_index_set, float bs_w, float ss_w);
	float Get_Avg_Dist(User* u, POI* p);
};

