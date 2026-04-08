 #pragma once
#include<vector>
#include "POI.h"
#include <map>
#include <set>
#include <bitset>
#include <unordered_map>
//static const int MAX_KEYS = 50;
class User
{
public:
	int id;
	
	
	std::map<int,float> inNeighbors;// follower id with thier influence 
	std::map<int, float> outNeighbors;// following id with thier influence 

	std::vector<double> Keys_Fsum;
	std::vector<double> Keys_Fmax;
	std::vector<int> Keys_visited_count;

	int ub_f_sum;
	float ub_f_avg;
	float ub_sub;
	float ub_w_in;
	float ub_w_out;
	int P_s_size;// Number of social network pivots
	std::vector<int> dist_P_s;// distance between user and all social network pivots;with size = P_s.size();  where index i points to spiv_i

	User(int _id,int _P_s_size, int total_n_keys);// need to specify P_r and P_s sizes
	~User();
	

	//For Pruning
	bool isPruned;


	std::map<int, int> checkin_locations;// point of interest that user visited with the frequesncy, always valid total valid not expired
	std::unordered_map<int, std::vector<int>> timed_checkins;

	int last_batch_updated = 0;
};
