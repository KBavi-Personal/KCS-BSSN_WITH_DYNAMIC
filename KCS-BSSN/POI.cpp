#include "POI.h"

POI::POI(int _id, int _cat_id, float coor_x, float coor_y, int num_keys, int _P_r_size,int totsl_n_keys) {
	id = _id;
	cat_id = _cat_id;
	coordinate_x = coor_x;
	coordinate_y = coor_y;
	num_of_keys = num_keys;
	keys = new int[num_keys];
	P_r_size = _P_r_size;
	const float INF = std::numeric_limits<float>::infinity();
	dist_P_r.resize(P_r_size, INF);

	Keys_Fsum.resize(totsl_n_keys,0.0);
	Keys_visited_count.resize(totsl_n_keys,0);


	isPruned = false;
}

POI::~POI() {
	if (keys != NULL && num_of_keys > 0) {
		delete[] keys;
		keys = NULL;
		num_of_keys = 0;
	}
	
}
bool POI::HasKey(int key)
{
	for (int k = 0; k < num_of_keys; k++)
		if (keys[k] == key)
			return true;
	return false;

}
bool POI::HasAnyKey(std::vector<int> Q)
{
	for (int key : Q)
		if (HasKey(key))return true;
	return false;
}




