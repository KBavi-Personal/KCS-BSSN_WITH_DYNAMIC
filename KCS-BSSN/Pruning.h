#pragma once
#include <vector>
#include "SampleData.h"
#include "TreeNode.h"
class User;
class Pruning
{
public:
	SampleData* sd;
	Pruning(SampleData *_sd);
	~Pruning();

	float Get_lb_avgdist(User* u, User* q, std::vector<int>& Q) const;

	//Prune a user
	bool Keyword_based_pruning(User *u, std::vector<int>& Q) const;
	bool Omega_based_pruning(User* u, int Omega) const;
	bool Pi_based_pruning(User* u, float Pi) const;
	bool Influnce_based_pruning(User* u, User* v, float theta) const;
	bool Structural_cohesiveness_pruning_sub(User* u, int k) const;
	bool Social_distance_based_pruning(User* u, User* q,int d) const;
	bool Spatial_distance_based_pruning(User* u, User* q, std::vector<int>& Q, float sigma) const;

	//Prune a Index Tree Node
	bool Keyword_based_pruning_index_node(TreeNode * node, std::vector<int>& Q, int d) const;
	bool Omega_based_pruning_index_node(TreeNode* node, std::vector<int>& Q, int Omega, int d) const;
	bool Pi_based_pruning_index_node(TreeNode* node, std::vector<int>& Q, float Pi, int d) const;
	bool Influnce_based_pruning_index_node(TreeNode* node, User* q, float theta) const;
	bool Structural_cohesiveness_pruning_sub_index_node(TreeNode* node, int k) const;
	bool Social_distance_based_pruning_index_node(TreeNode* node, User* q, int d) const;

private:

	float LB_poi_poi(int poiA, int poiB) const;
};

