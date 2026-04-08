#pragma once
#include "SampleData.h"
#include <unordered_set>
#include <unordered_map>
#include "Indexing.h"
#include "CommStru.h"
#include "Helpers.h"
class Tree;
class Ref2;
class DynamicApproch
{
public:
	DynamicApproch();
	~DynamicApproch();

	//dynamic update related TO TREE
	SampleData* sd;
	Indexing* inx;
	std::vector<TreeNode*> all_nodes; 
	std::vector<TreeNode*> user_leaf;
	void UpdateTreeAfterBatch_WithOut_Migration(const std::vector<int>& updated_users);
	void RecomputeLeafMetrics( TreeNode* leaf);
	void RecomputeInternalMetrics( TreeNode* node);

	//Dynamic update with migration	
	std::vector<User*> P_index_vec;
	double ComputeNodePivotQuality(TreeNode* node, User* piv, float bs_w, float ss_w);
	std::vector<TreeNode*> GetAllLeafNodes() const;
	std::vector<TreeNode*> GetNodesAtDepth(int depth) const;
	bool MigrateUserToLeafNode( int user_id, float bs_w, float ss_w, float rs_w, double stability_margin, const std::vector<TreeNode*>& leaves);
	bool MigrateNodeToParentNode( TreeNode* node, float bs_w, float ss_w, double stability_margin);
	void UpdateTreeAfterBatch_With_Migration(float bs_w, float ss_w, float rs_w, double stability_margin, const std::vector<int>& updated_users);
	std::unordered_map<int, std::vector<TreeNode*>> nodes_by_depth;
	std::vector<TreeNode*> leaf_nodes_cached;
	std::map<int, std::unordered_set<TreeNode*>> touched_nodes;
	void RebuildUserWindowedCheckins(User* u, int current_day, int window_days);
	void RebuildUserKeyStats(User* u);
	void RecomputeUserUb(User* u);


	std::vector<UserVisit> Picked_users_get_visits(std::vector<int> picked_users, int max_visits, int min_day, int max_day);
	bool Validate_Communities_After_Repair(std::string filename, Ref2* ref, int q_index,
		std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value);
	void Add_affected_Users_With_POI_to_Remaining_one(int u_id, CommStru& comm);
	void Add_affected_Users_With_POI_to_Remaining_all(std::vector<int> affected_users, CommStru& comm);
	void Add_random_timestamps_for_first_time(std::string filename);
	void GenerateOfflineData1(int pivot_iterations, std::string filename, int distribution_type);
	std::vector<CommStru> Load_Communities(std::string filename);
	void Reset_Remaining_By_Community(CommStru& comm);
	void Picked_users_Add_visit(const std::vector<UserVisit>& user_visits);
	std::vector<int> Picked_users_Filter_visit( std::vector<int> picked_users, int max_day, int window_days);
	std::string FindCommuntiesPrintToFile_Defults( Tree* tree, std::string filename, Ref2* ref,
		std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value);
		std::string Run_Dynamic_Temporal_Delete_visits_by_user_number(SampleData* sd, int user_no_updates, std::vector<int> picked_users, int max_day, int window_days, std::string filename,
		std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value,float bs_w, float ss_w, float rs_w, double stability_margin, std::unordered_set<int> comm_indeies);
	std::string Run_Dynamic_Temporal_Insert_visits(SampleData* sd, int user_no_updates, std::vector<int> picked_users, int max_visits, int min_day, int max_day, std::string filename,
		std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, 
		float theta_value, float sigma_value, float bs_w, float ss_w, float rs_w, double stability_margin, std::unordered_set<int> comm_indeies);
	void Build_tree_and_save_withComm( std::string filename, bool loadTree,
		std::vector<int> Q, int k_value, int d_value, float omega_value, float pi_value, float theta_value, float sigma_value);
	std::unordered_set<int> Select_random_comm_index_to_test(std::string filename, int num_to_select);




};

