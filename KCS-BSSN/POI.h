#pragma once
#include <string>
#include <vector>
class POI
{
public:
	int id;
	int cat_id;
	float coordinate_x;// location of the POI
	float coordinate_y;// location of the POI
	int num_of_keys;
	int* keys; //index of the key; map to keys from keys description list
	int my_closest_intesect_point;

	int P_r_size;// Number of road network pivots
	std::vector<float> dist_P_r;// distance between user and all spatial network pivots;with size = P_r.size();  where index i points to rpiv_i
	POI(int _id, int cat_id, float coor_x, float coor_y, int num_keys, int P_r_size, int total_n_keys);
	~POI();
	bool HasKey(int key);
	bool HasAnyKey(std::vector<int> Q);
	std::vector<double> Keys_Fsum;
	std::vector<int> Keys_visited_count;

	//for pruning 
	bool isPruned;
};

