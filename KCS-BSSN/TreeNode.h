#pragma once
#include <map>
#include <vector>
#include <string>

class User;
class SampleData;
class Tree;
class TreeNode
{
public:
	int id;
	int depth = 0;
	std::vector<User*> users;
	Tree* my_tree;
	bool is_leaf;
	std::vector<TreeNode*> children;//tree nodes
	bool isroot;
	int total_n_keys;

	std::vector<double> Keys_ub_Fsum;
	std::vector<double> Keys_ub_Fmax;
	//std::map<int, int> Key_f_sum_map;//  keyword, f_sum
	float ub_sub;
	float ub_w_in;
	float ub_w_out;
	int P_r_size;// Number of road network pivots
	int P_s_size;// Number of social network pivots
	float* min_max_dist_P_s;// distance between user and all social network pivots;with size=2*P_s; where item i and i+1 points to min and max avgdist to spiv_i

	TreeNode(int level, int P_s_size, int P_r_size, Tree *tree,bool isLeaf, int pivot_user_index, int _total_n_keys);// need to specify P_s sizes;
	~TreeNode();

	int myUserIndexForPivot;
	TreeNode* parent = nullptr;

	double min_bs = 0.0, max_bs = 0.0;
	double min_rs = 0.0, max_rs = 0.0;
	double min_ss = 0.0, max_ss = 0.0;
	

};

