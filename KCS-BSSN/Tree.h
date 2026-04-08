#pragma once
#include "TreeNode.h"
#include "Indexing.h"
#include <unordered_set>
#include <unordered_map>


class TreeNode;
class User;
class SampleData;
class Indexing;
class Tree
{
public:
	SampleData* sd;
	Indexing* inx;
	TreeNode* root;
	int num_of_nodes;
	std::vector<TreeNode*> Nodes;
	std::vector<TreeNode*> all_nodes;
	std::vector<User*> P_index_vec_for_nodes;
	void BuildTree(int iteration, int  tree_level_division, int piv_size, float bs_w, float ss_w, float rs_w);
	void ButtomUp(int tree_level_division, int piv_selction_iter, float bs_w, float ss_w, float rs_w);;

	void BulidLeafNodes();
	TreeNode* Create_Leaf_Node(Indexing* inx, int pivot_user_index, const std::vector<User*>& group);
	TreeNode* Create_Non_Leaf_Node(int level, std::vector<TreeNode*> Nodes, int pivot_user_index);

	Tree(SampleData* sd);
	~Tree();

	bool KCS_BSSN_Query_Answer(std::vector<int> Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index, Sextuple<bool, bool, bool, bool, bool, bool, bool> m);
	std::vector<TreeNode*> user_leaf;
	
};

