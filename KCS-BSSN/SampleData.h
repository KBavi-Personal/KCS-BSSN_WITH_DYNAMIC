#pragma once
#include "BipartiteNetwork.h"
#include <set>
#include <map>
#include <fstream>

class SampleData
{
public:
	BipartiteNetwork* bn;
	int P_r_size;//number of pivots on road netwrok
	int P_s_size;//number of pivots on social netwrok
	int no_of_POI;
	int max_x;
	int max_y;
	std::vector<User*> P_s_vec;
	std::vector<int> P_r_vec;//index of intersection point

	std::vector<POI*> poi_vec;
	std::vector<float> intersection_points;
	int total_n_of_intersection_points;
	int toatl_n_of_used_keys;
	std::map<int, std::map<int, float>> roadNetworkAdjMtrx;//i -> j, dist


	//float* intersection_points_roaddist;// shortest road dist between i*total_n_of_intersection_points+j
	int larget_f_sum;
	float larget_f_avg;

	int no_remining_users;// Number of remining users 
	int* remining_users;// Reminng users ids
	int no_remining_POI;// Number of remining POIs 
	int* remining_POI;// Reminng POI ids

	SampleData(int num_of_users, int total_n_of_POI, int no_p_r, int no_p_s, int max_x, int max_y, int total_n_of_intersection_points, int _toatl_n_of_used_keys);
	~SampleData();
	void Generate_POI(int max_no_of_keys_per_p, int distribution_type);
	void Generate_intersection_points(bool isUniform);
	void Generate_edges_road_Network(bool write_to_file);
	void Generate_Social_Network_Random_Users(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency);
	void Generate_Social_Network_Random_edges(int _max_no_follower, bool IsUniform);



	long double Generate_road_network(bool isUniform);
	float Get_avg_dist_user_poi(User* u, int p_index);
	float dist_r(int p1, int p2);

	void Print_poi(bool IsWrite, std::string filename);
	void Print_SN(bool IsWrite, std::string filename);
	void Print_intersection_points(bool IsWrite);
	void Print_road_network_nodes(bool IsWrite, int no);
	void Print_road_network_edges(bool IsWrite, int no);
	void Print_dist_road_network(bool IsWrite);


	//real dataset
	double Load_Twitter_Social_Network_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency);
	double Load_Epinions_Social_Network_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency);

	void Load_cal_RN_POI_BN();
	void Load_Synthatic_RN();
	void Load_SN_BN(std::string filename);
	void Load_POI(std::string filename);
	void Load_BN(std::string dataset, bool IsRealDataset, int distribution_type);// uniform=1, guassian =2 , zipf=3


	void Print_To_File(std::string filename, std::string text);
	void Print_Graph_Edges_ToFile();
	void Reset_POI_pruned_flag()
	{
		for (auto p : poi_vec)
			p->isPruned = false;
	}
	void Reset_Users_pruned_flag()
	{
		for (int i = 0; i < bn->no_of_users; i++)
			bn->users_vec[i]->isPruned = false;
	}
	void Reset_Remain_users_flag()
	{
		for (int i = 0; i < no_remining_users; i++)
			bn->users_vec[remining_users[i]]->isPruned = false;
	}
	void Generate_Social_Network_Random_Users_add_extra(int new_user_size);
	void Assign_poi_to_users(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency, int new_user_size, int  start_index);
	void fix_dblp_ungraph();
	double Load_dblp_ungraph_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency);

	void generate_Socail_Netword_with_k_d_truss(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequency, bool isUniformFrequency, std::string filename, std::vector<int> users_size_vec, int _min_no_follower, int _max_no_follower, int distribution_type);  // uniform=1, guassian =2 , zipf=3

	void Print_SN_timestamp(bool IsWrite, std::string filename);



};