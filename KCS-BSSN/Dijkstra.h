#pragma once
#include<vector>
#include<set>

class BipartiteNetwork;
class SampleData;
class User;
class Dijkstra {
public:
	SampleData* sd;
	Dijkstra(SampleData* sd);
	~Dijkstra();

	float Calculate_shotest_dist_road(SampleData* sd, int src_location_index, int destination_location_index);
	void Calculate_shotest_dist_road_All_intersection_points(SampleData* sd);
	std::vector<float> MaxInfluenceFromSource_All_OutNeighbors(int src_user_Id, float theta);

	void Calculate_shotest_dist_road_s_to_all(int pr_index);
};