#include "SampleData.h"
#include "BipartiteNetwork.h"
#include <iostream>
#include <cmath>
#include "BFS.h"
#include "MinHeap.h"
#include "Dijkstra.h"
#include "rand.h"
#include <unordered_set>
#include <random>
#include "SocialNetworkGraphKTruss.h"


SampleData::SampleData(int num_of_users, int total_n_of_POI, int no_p_r, int no_p_s, int _max_x, int _max_y, int _total_n_of_intersection_points, int _toatl_n_of_used_keys)
{
	P_r_size = no_p_r;
	P_s_size = no_p_s;
	bn = new BipartiteNetwork();
	bn->no_of_users = num_of_users;
	no_of_POI = total_n_of_POI;
	max_x = _max_x;
	max_y = _max_y;
	total_n_of_intersection_points = _total_n_of_intersection_points;
	toatl_n_of_used_keys = _toatl_n_of_used_keys;

	intersection_points.resize(2 * total_n_of_intersection_points);

	no_remining_users = 0;
	remining_users = new int[num_of_users];
	no_remining_POI = 0;
	remining_POI = new int[total_n_of_POI];

	larget_f_sum = 0;
	larget_f_avg = 0.0f;
}
SampleData::~SampleData() {
	delete[] remining_users;
	delete[] remining_POI;
	delete bn;
	for (POI* p : poi_vec)
		delete p;

}
float cal_distance(float x1, float x2, float y1, float y2) {
	float dist = sqrt(pow(abs(x2 - x1), 2) + pow(abs(y2 - y1), 2));
	return dist;
}
void SampleData::Generate_intersection_points(bool isUniform) {
	float p_x, p_y;
	for (int i = 0; i < total_n_of_intersection_points; i++) {
		p_x = isUniform ? random1::uniform_int_pos(0, max_x * 10) / 10.0 : random1::zipf_pos(0, max_x * 10, 0.8) / 10.0;
		p_y = isUniform ? random1::uniform_int_pos(0, max_y * 10) / 10.0 : random1::zipf_pos(0, max_y * 10, 0.8) / 10.0;

		intersection_points.at(2 * i) = p_x;
		intersection_points.at(2 * i + 1) = p_y;

	}
}
void SampleData::Generate_edges_road_Network(bool write_to_file) {
	//generate  PLANNER graph using Gabriel graph algorithm
	for (int i = 0; i < total_n_of_intersection_points; i++)
		for (int j = 0; j < total_n_of_intersection_points; j++) {
			if (i == j) continue;
			float mid_x = (intersection_points.at(2 * i) + intersection_points.at(2 * j)) / 2;
			float mid_y = (intersection_points.at(2 * i + 1) + intersection_points.at(2 * j + 1)) / 2;
			float radius = cal_distance(mid_x, intersection_points.at(2 * i), mid_y, intersection_points.at(2 * i + 1));
			//if (uniform(0, 10) < 6)continue;// add probabity to connect 0.6
			bool connect = true;
			for (int k = 0; k < total_n_of_intersection_points; k++) {
				if (k == i || k == j) continue;
				float distance = cal_distance(mid_x, intersection_points.at(2 * k), mid_y, intersection_points.at(2 * k + 1));
				if (distance <= radius) connect = false;
			}
			if (connect) {
				float rd_dist = cal_distance(intersection_points.at(2 * i), intersection_points.at(2 * j), intersection_points.at(2 * i + 1), intersection_points.at(2 * j + 1));
				roadNetworkAdjMtrx[i].insert(std::make_pair(j, rd_dist));

			}
		}
	if (write_to_file)
		Print_road_network_edges(true, total_n_of_intersection_points);
}
long double SampleData::Generate_road_network(bool isUniform)
{
	long double index_time;
	long double ind_st, ind_ed;
	ind_st = clock();
	Generate_intersection_points(isUniform);
	ind_ed = clock();
	//index 0, Generate_intersection_points_time
	printf("Generate_intersection_points time= %f\n", (ind_ed - ind_st) / 1000.0);
	index_time = (ind_ed - ind_st) / 1000.0;
	Print_road_network_nodes(true, total_n_of_intersection_points);

	ind_st = clock();
	Generate_edges_road_Network(false);
	ind_ed = clock();
	//index 1, Generate_edges_road_Network_time
	index_time += (ind_ed - ind_st) / 1000.0;
	printf("Generate_edges_road_Network time= %f\n", (ind_ed - ind_st) / 1000.0);

	return index_time;
}
int get_random_key_index(int distribution_type, int max, int min) {
	int random_key_index;
	if (distribution_type == 3) {  // zipf
		random_key_index = random1::zipf_pos(min, max, 0.8);
	}
	else if (distribution_type == 2) { // gaussian
		random_key_index = random1::gaussian_pos(min, max, max / 2.0, max / 5.0);
	}
	else { // default to uniform
		random_key_index = random1::uniform_int_pos(min, max);
	}
	return random_key_index;
}
void SampleData::Generate_POI(int max_no_of_keys_per_p, int distribution_type) { // uniform=1, guassian =2 , zipf=3
	Load_Synthatic_RN();
	float p_x, p_y;
	int p_no_of_keys, random_key_index;
	for (int i = 0; i < no_of_POI; i++) {
		p_x = random1::uniform_int_pos(0, max_x * 10) / 10.0;
		p_y = random1::uniform_int_pos(0, max_y * 10) / 10.0;
		p_no_of_keys = random1::uniform_int_pos(3, max_no_of_keys_per_p);
		POI* p = new POI(i, -1, p_x, p_y, p_no_of_keys, P_r_size, toatl_n_of_used_keys);
		std::set<int> temp_keys;
		for (int j = 0; j < p_no_of_keys; j++) {


			random_key_index = get_random_key_index(distribution_type, toatl_n_of_used_keys - 1, 0);
			while (temp_keys.count(random_key_index) == 1)
				random_key_index = get_random_key_index(distribution_type, toatl_n_of_used_keys - 1, 0);
			temp_keys.insert(random_key_index);
			p->keys[j] = random_key_index;
		}
		float min_dist = std::numeric_limits<float>::max();
		for (int k = 0; k < total_n_of_intersection_points; k++) {
			float dist = cal_distance(p->coordinate_x, intersection_points.at(2 * k), p->coordinate_y, intersection_points.at(2 * k + 1));

			if (min_dist > dist) {
				min_dist = dist;
				p->my_closest_intesect_point = k;
			}
		}
		poi_vec.push_back(p);
	}
	poi_vec.shrink_to_fit();
}
void SampleData::Generate_Social_Network_Random_Users(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency)
{
	int user_no_of_checkin_loc, poi_index;

	for (int i = 0; i < bn->no_of_users; i++) {
		User* newUser = new User(i, P_s_size, toatl_n_of_used_keys);
		user_no_of_checkin_loc = random1::uniform_int_pos(1, max_no_checkin_loc);

		bool availble_euclidean_dist_btwn_user_POI = true;
		int first_poi_index = -1;
		for (int j = 0; j < user_no_of_checkin_loc; j++)
		{
			poi_index = random1::uniform_int_pos(0, poi_vec.size() - 1);
			if (first_poi_index != -1) {
				int counter = 0;
				bool isinList = false;
				for (auto ch_loc : newUser->checkin_locations) {
					if (poi_index == ch_loc.first) {
						isinList = true;
						break;
					}
				}

				while (isinList || cal_distance(poi_vec[first_poi_index]->coordinate_x, poi_vec[poi_index]->coordinate_x, poi_vec[first_poi_index]->coordinate_y, poi_vec[poi_index]->coordinate_y) > max_euclidean_dist_btwn_user_POI)
				{
					poi_index = random1::uniform_int_pos(0, poi_vec.size() - 1);
					counter++;
					isinList = false;
					for (auto ch_loc : newUser->checkin_locations) {
						if (poi_index == ch_loc.first) {
							isinList = true;
							break;
						}

					}
					//if it seached for 1000 times and didn't find any neighbors then pass
					if (counter > 1000) {
						availble_euclidean_dist_btwn_user_POI = false;
						break;
					}
				}
				if (!availble_euclidean_dist_btwn_user_POI)
					break;
			}
			else {
				first_poi_index = poi_index;
			}
			int frequncy = 0;
			while (frequncy <= 0)
				frequncy = isUniformFrequency ? random1::uniform_int_pos(1, max_frequncy) : random1::gaussian_pos(1, max_frequncy, max_frequncy / 2.0, max_frequncy / 5.0);// assign random frequency to checkin_loc
			newUser->checkin_locations[poi_index] = frequncy;
		}
		bn->users_vec.push_back(newUser);
	}
	bn->users_vec.shrink_to_fit();
}
void SampleData::Generate_Social_Network_Random_edges(int _max_no_follower, bool IsUniform)
{
	for (int i = 0; i < bn->no_of_users; i++) {
		int no_of_Neighbors = random1::uniform_int_pos(1, _max_no_follower);

		no_of_Neighbors = no_of_Neighbors > 0 ? no_of_Neighbors : 1;
		for (int j = 0; j < no_of_Neighbors; j++) {
			int random_user_id = -1;
			while (random_user_id<0 || random_user_id>bn->no_of_users - 1)
				random_user_id = IsUniform ? random1::uniform_int_pos(0, bn->no_of_users - 1) : random1::gaussian_pos(0, bn->no_of_users, (bn->no_of_users - 1) / 2.0, (bn->no_of_users - 1) / 5.0);

			if (random_user_id == i) { i--; continue; }// if it is same as user id, ignore and redo
			float random_influnce = random1::gaussian_pos(1, 10, 10 / 2.0, 10) / 10;
			bn->users_vec[i]->outNeighbors[random_user_id] = random_influnce;
			bn->users_vec[random_user_id]->inNeighbors[i] = random_influnce;
			bn->Set_edge_sup(i, random_user_id, 0);

		}
	}
}
float SampleData::Get_avg_dist_user_poi(User* u, int p_index)
{
	float avg_dist = 0;
	for (std::pair<int, int> check_in_loc : u->checkin_locations)
	{
		avg_dist += dist_r(check_in_loc.first, p_index);
	}
	return avg_dist / u->checkin_locations.size();
}
float SampleData::dist_r(int p1, int p2)
{
	float dist;
	Dijkstra* dijkstra = new Dijkstra(this);
	dist = dijkstra->Calculate_shotest_dist_road(this, p1, p2);
	delete dijkstra;
	return dist;
}
void SampleData::Print_intersection_points(bool IsWrite)
{
	printf("\Print_intersection_points\n");
	if (IsWrite) {
		std::ofstream myfile;
		myfile.open("intersectionPoints.txt");
		for (int i = 0; i < total_n_of_intersection_points; i++)
			myfile << i << " " << intersection_points.at(2 * i) << " " << intersection_points.at(2 * i + 1) << "\n";
		myfile.close();
	}
	else {
		int i = 0;
		for (int i = 0; i < total_n_of_intersection_points; i++)
			printf("(id=%d - loc(%f,%f) \n", i, intersection_points.at(2 * i), intersection_points.at(2 * i + 1));
	}
}
void SampleData::Print_poi(bool IsWrite, std::string filename)
{
	printf("\POI\n");
	if (!IsWrite) {

		for (POI* p : poi_vec) {
			printf("(id=%d - loc(%f,%f) closest_intersect=%d (%f,%f) %d-keys=(", p->id, p->coordinate_x, p->coordinate_y,
				p->my_closest_intesect_point, intersection_points.at(2 * p->my_closest_intesect_point),
				intersection_points.at(2 * p->my_closest_intesect_point + 1), p->num_of_keys);

			printf(")\n");
		}
	}
	else {
		std::string file_name = filename;
		std::ofstream myfile;
		myfile.open(file_name);
		for (POI* p : poi_vec) {
			myfile << p->id << " " << p->cat_id << " " << p->coordinate_x << " " << p->coordinate_y << " " <<
				p->my_closest_intesect_point << " " << p->num_of_keys;
			for (int i = 0; i < p->num_of_keys; i++)
			{
				myfile << " " << p->keys[i];
			}
			myfile << "\n";
		}
		myfile.close();
	}
}
void SampleData::Print_SN(bool IsWrite, std::string filename)
{
	printf("\SN\n");
	if (!IsWrite) {
		filename = filename + "_sn.txt";
		std::string text = bn->print();
		std::ofstream myfile;
		myfile.open(filename);
		myfile << text;
		myfile.close();
	}
	else {
		std::string file_name = filename + "ip_" + std::to_string(total_n_of_intersection_points) + "_users_" + std::to_string(bn->no_of_users) + "_poi_" + std::to_string(no_of_POI) + "_sn.txt";
		std::ofstream myfile;
		myfile.open(file_name);
		for (int i = 0; i < bn->no_of_users; i++) {
			User* u = bn->users_vec[i];
			myfile << u->id;
			for (auto& u_ch_in : u->checkin_locations) {
				myfile << " " << u_ch_in.first << ":" << u_ch_in.second;
			}
			myfile << "&";
			for (auto& u_in_n : u->inNeighbors) {
				myfile << " " << u_in_n.first << ":" << u_in_n.second;
			}
			myfile << "#";
			for (auto& u_out_n : u->outNeighbors) {
				myfile << " " << u_out_n.first << ":" << u_out_n.second;
			}
			myfile << "\n";
		}
		myfile.close();
	}
}
void SampleData::Print_road_network_nodes(bool IsWrite, int no)
{
	printf("\Road Network\n");
	std::string file_name = std::to_string(no) + "_nodes.txt";
	if (IsWrite) {
		std::ofstream myfile;
		myfile.open(file_name);
		for (int i = 0; i < total_n_of_intersection_points; i++)
			myfile << i << " " << intersection_points.at(2 * i) << " " << intersection_points.at(2 * i + 1) << "\n";
		myfile.close();
	}

}
void SampleData::Print_road_network_edges(bool IsWrite, int no)
{
	printf("\Road Network\n");
	std::string file_name = std::to_string(no) + "_edges.txt";
	if (IsWrite) {
		std::ofstream myfile;
		myfile.open(file_name);
		for (auto r : roadNetworkAdjMtrx)
			for (auto s : r.second) {
				myfile << r.first << "-" << s.first << "=" << s.second << "\n";

			}
		myfile.close();
	}
	else {
		for (auto r : roadNetworkAdjMtrx)
			for (auto s : r.second) {
				printf("%d %d %f\n", r.first, s.first, s.second);
			}
		printf("\n");
	}
}
void SampleData::Print_dist_road_network(bool IsWrite)
{
	std::ofstream myfile;
	if (IsWrite)
		myfile.open("distRoadNetwork.txt");


	if (!IsWrite)
		printf("dist Road Network map\n");

	if (!IsWrite)
		printf("\dist Road Network\n");
	for (int i = 0; i < total_n_of_intersection_points; i++)
		for (int j = 0; j < total_n_of_intersection_points; j++)
		{
			if (IsWrite)
				myfile << i << "-" << j << "=" << 0 << "\n";
			else
				printf("%d-%d\n", i, j);
		}
	myfile.close();
}
void fix_twitter_sn_network() {
	std::ifstream fileNode("dataset/sn/twitter.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open twitter file!\n";
		return;
	}
	std::string delimiter = " ";
	std::string line;


	std::set<int> nodes;
	while (getline(fileNode, line))
	{
		int id = std::stoi(line.substr(0, line.find(delimiter)));
		if (nodes.count(id) == 0)
			nodes.insert(id);
		line.erase(0, line.find(delimiter) + 1);
		id = std::stoi(line.substr(0, line.find(delimiter)));
		if (nodes.count(id) == 0)
			nodes.insert(id);
	}
	fileNode.close();
	std::ofstream myfile1;
	myfile1.open("twitter_mapping.txt");
	std::map<int, int> maap;
	int i = 0;
	for (int node : nodes) {
		maap[node] = i;
		myfile1 << node << ":" << i << "\n";
		i++;

	}
	myfile1.close();

	std::ifstream fileNode2("dataset/sn/twitter.txt");
	std::ofstream myfile;
	myfile.open("twitter_by_index.txt");
	while (getline(fileNode2, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));

		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		myfile << maap[id1] << " " << maap[id2] << "\n";
	}
	fileNode2.close();
	myfile.close();
}
double SampleData::Load_Twitter_Social_Network_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency)
{
	Generate_Social_Network_Random_Users(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequncy, isUniformFrequency);

	std::ifstream fileNode("dataset/sn/twitter_by_index.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file!\n";
		return 0;
	}
	std::string delimiter = " ";
	std::string line;
	while (getline(fileNode, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		float random_influnce = random1::gaussian_pos(1, 10, 10 / 2.0, 10) / 10;
		bn->users_vec[id1]->outNeighbors[id2] = random_influnce;
		bn->users_vec[id2]->inNeighbors[id1] = random_influnce;
		bn->Set_edge_sup(id1, id2, 0);
	}
	fileNode.close();
	//Print_SN(true, "twitter_sn");
	return 0;
}
double SampleData::Load_Epinions_Social_Network_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency)
{
	Generate_Social_Network_Random_Users(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequncy, isUniformFrequency);

	std::ifstream fileNode("dataset/sn/soc-Epinions1.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file!\n";
		return 0;
	}
	std::string delimiter = "\t";
	std::string line;
	while (getline(fileNode, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		float random_influnce = random1::gaussian_pos(1, 10, 10 / 2.0, 10) / 10;
		bn->users_vec[id1]->outNeighbors[id2] = random_influnce;
		bn->users_vec[id2]->inNeighbors[id1] = random_influnce;
		bn->Set_edge_sup(id1, id2, 0);

	}
	fileNode.close();

	return 0;
}
void SampleData::Load_cal_RN_POI_BN() {

	std::ifstream fileNode("dataset/cal/calcnode.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file calcnode!\n";
		return;
	}
	std::string delimiter = " ";
	std::string line;
	int i = 0;
	while (getline(fileNode, line))
	{
		line.erase(0, line.find(delimiter) + 1);
		intersection_points.at(2 * i) = std::stof(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		intersection_points.at(2 * i + 1) = std::stof(line.substr(0, line.find(delimiter)));
		i++;
	}
	fileNode.close();


	std::ifstream fileEdge("dataset/cal/calcedge.txt");

	if (!fileEdge.is_open())
	{
		std::cout << "Could not open file calcedge!\n";
		return;
	}
	int st_node, end_node;
	float dist;
	while (getline(fileEdge, line))
	{
		line.erase(0, line.find(delimiter) + 1);
		st_node = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		end_node = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		dist = std::stof(line.substr(0, line.find(delimiter)));
		roadNetworkAdjMtrx[st_node].insert(std::make_pair(end_node, dist));
		roadNetworkAdjMtrx[end_node].insert(std::make_pair(st_node, dist));
	}
	fileEdge.close();
	Load_POI("dataset/cal/calpoi_info.txt");
}
void SampleData::Load_POI(std::string filename) {
	std::string delimiter = " ";
	std::string line;

	std::ifstream filePOI(filename);

	if (!filePOI.is_open())
	{
		std::cout << "Could not open file " << filename << " !\n";
		return;
	}
	float p_x, p_y;
	int p_id, cat_id, p_no_of_keys, my_closest_intesect_point;
	//format id cat_id coordinate_x coordinate_y my_closest_intesect_point num_of_keys keys[0] keys[1] keys[2] ...
	while (getline(filePOI, line))
	{
		p_id = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		cat_id = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		p_x = std::stof(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		p_y = std::stof(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		my_closest_intesect_point = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		p_no_of_keys = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);


		POI* p = new POI(p_id, cat_id, p_x, p_y, p_no_of_keys, P_r_size, toatl_n_of_used_keys);
		p->my_closest_intesect_point = my_closest_intesect_point;
		for (int j = 0; j < p_no_of_keys; j++) {
			p->keys[j] = std::stoi(line.substr(0, line.find(delimiter)));
			line.erase(0, line.find(delimiter) + 1);
		}
		poi_vec.push_back(p);
	}
	filePOI.close();
	poi_vec.shrink_to_fit();
}
void SampleData::Load_Synthatic_RN()
{

	std::string file_name = std::to_string(total_n_of_intersection_points) + "_nodes.txt";
	std::ifstream fileRnodes("dataset/synthatic/rn/" + file_name);
	if (!fileRnodes.is_open())
	{
		std::cout << "Could not open file! " << file_name << "\n";
	}
	std::string line;
	int i;
	while (getline(fileRnodes, line))
	{
		i = std::stoi(line.substr(0, line.find(" ")));
		line.erase(0, line.find(" ") + 1);
		intersection_points.at(2 * i) = std::stof(line.substr(0, line.find(" ")));
		line.erase(0, line.find(" ") + 1);
		intersection_points.at(2 * i + 1) = std::stof(line.substr(0, line.find(" ")));
	}
	fileRnodes.close();


	file_name = std::to_string(total_n_of_intersection_points) + "_edges.txt";
	std::ifstream fileRedgest("dataset/synthatic/rn/" + file_name);
	if (!fileRedgest.is_open())
	{
		std::cout << "Could not open file! " << file_name << "\n";
	}

	int st_node, end_node;
	float dist;
	while (getline(fileRedgest, line))
	{
		st_node = std::stoi(line.substr(0, line.find("-")));
		line.erase(0, line.find("-") + 1);
		end_node = std::stoi(line.substr(0, line.find("=")));
		line.erase(0, line.find("=") + 1);
		dist = std::stof(line.substr(0, line.find(" ")));
		roadNetworkAdjMtrx[st_node].insert(std::make_pair(end_node, dist));
		roadNetworkAdjMtrx[end_node].insert(std::make_pair(st_node, dist));
	}
	fileRedgest.close();
}
void SampleData::Load_SN_BN(std::string dataset)
{
	std::string file_name = dataset + "ip_" + std::to_string(total_n_of_intersection_points) + "_users_" + std::to_string(bn->no_of_users) + "_poi_" + std::to_string(no_of_POI) + "_sn.txt";
	std::ifstream fileNode("dataset/bn/" + file_name);

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file! " << file_name << "\n";
	}
	std::string line, subline, subsubline;
	//format u->id [u_ch_in:frequency]& [inNeighbor:influence]# [outNeighbors::influence]
	while (getline(fileNode, line))
	{
		int u_id = std::stoi(line.substr(0, line.find(" ")));
		line.erase(0, line.find(" ") + 1);

		User* newUser = new User(u_id, P_s_size, toatl_n_of_used_keys);
		bn->users_vec.push_back(newUser);

		//checkin loc
		subline = line.substr(0, line.find("&"));
		line.erase(0, line.find("&") + 1);
		//subline.erase(0, subline.find(" ") + 1);
		while (subline.length() > 1) {
			subsubline = subline.substr(0, subline.find(" "));
			subline.find(" ") != -1 ? subline.erase(0, subline.find(" ") + 1) : subline = "";

			int id = std::stoi(subsubline.substr(0, subsubline.find(":")));
			subsubline.erase(0, subsubline.find(":") + 1);
			int value = std::stoi(subsubline);
			bn->users_vec[u_id]->checkin_locations[id] = value;
		}

		//InNeighbors
		subline = line.substr(0, line.find("#"));
		line.erase(0, line.find("#") + 1);
		subline.erase(0, subline.find(" ") + 1);
		while (subline.length() > 1) {
			subsubline = subline.substr(0, subline.find(" "));
			subline.find(" ") != -1 ? subline.erase(0, subline.find(" ") + 1) : subline = "";

			int id = std::stoi(subsubline.substr(0, subsubline.find(":")));
			subsubline.erase(0, subsubline.find(":") + 1);
			bn->users_vec[u_id]->inNeighbors[id] = std::stof(subsubline);
		}
		//OutNeighbors
		line.erase(0, line.find(" ") + 1);
		subline = line;
		while (subline.length() > 1) {
			subsubline = subline.substr(0, subline.find(" "));
			subline.find(" ") != -1 ? subline.erase(0, subline.find(" ") + 1) : subline = "";

			int id = std::stoi(subsubline.substr(0, subsubline.find(":")));
			subsubline.erase(0, subsubline.find(":") + 1);
			bn->users_vec[u_id]->outNeighbors[id] = std::stof(subsubline);
			bn->Set_edge_sup(u_id, id, 0);
		}
	}
	fileNode.close();
	bn->users_vec.shrink_to_fit();

}
void SampleData::Load_BN(std::string dataset, bool IsRealDataset, int distribution_type)// uniform=1, guassian =2 , zipf=3
{
	if (IsRealDataset) {
		Load_cal_RN_POI_BN();
	}
	else {
		Load_Synthatic_RN();
		std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";
		Load_POI("dataset/poi/" + distribution_type_name + "ip_" + std::to_string(total_n_of_intersection_points) + "_poi_" + std::to_string(no_of_POI) + ".txt");
	}
	Load_SN_BN(dataset);
}
void SampleData::Print_To_File(std::string filename, std::string text)
{
	std::string file_name = filename;

	std::ofstream myfile;
	myfile.open(file_name, std::ios_base::app);
	myfile << text << "\n";
	myfile.close();

}
void SampleData::Print_Graph_Edges_ToFile()
{
	std::string file_name = "SN_graph_edges_" + std::to_string(bn->no_of_users) + ".txt";

	std::ofstream myfile;
	myfile.open(file_name);
	for (auto u : bn->users_vec) {
		for (auto e : u->outNeighbors) {
			myfile << u->id << " " << e.first << " " << e.second << "\n";
		}
	}
	myfile.close();
}
void SampleData::Generate_Social_Network_Random_Users_add_extra(int new_user_size)
{
	int start_index = bn->users_vec.size();
	bn->no_of_users = new_user_size;
	for (int i = start_index; i < bn->no_of_users; i++) {
		User* newUser = new User(i, P_s_size, toatl_n_of_used_keys);
		bn->users_vec.push_back(newUser);
	}
	bn->users_vec.shrink_to_fit();
}
void SampleData::Assign_poi_to_users(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency, int new_user_size, int start_index)
{

	int user_no_of_checkin_loc, poi_index;
	for (int i = start_index; i < bn->no_of_users; i++) {
		user_no_of_checkin_loc = random1::uniform_int_pos(1, max_no_checkin_loc);

		bool availble_euclidean_dist_btwn_user_POI = true;
		int first_poi_index = -1;
		for (int j = 0; j < user_no_of_checkin_loc; j++)
		{
			poi_index = random1::uniform_int_pos(0, poi_vec.size() - 1);
			//while (first_poi_index != -1 && intersection_points_roaddist[poi_vec[first_poi_index]->my_closest_intesect_point * total_n_of_intersection_points + poi_vec[poi_index]->my_closest_intesect_point] > max_euclidean_dist_btwn_user_POI)
			if (first_poi_index != -1) {
				int counter = 0;
				bool isinList = false;
				for (auto ch_loc : bn->users_vec[i]->checkin_locations) {
					if (poi_index == ch_loc.first) {
						isinList = true;
						break;
					}
				}
				// assign checkin loc within max_euclidean_dist_btwn_user_POI radius

				while (isinList || cal_distance(poi_vec[first_poi_index]->coordinate_x, poi_vec[poi_index]->coordinate_x, poi_vec[first_poi_index]->coordinate_y, poi_vec[poi_index]->coordinate_y) > max_euclidean_dist_btwn_user_POI)
				{
					poi_index = random1::uniform_int_pos(0, poi_vec.size() - 1);
					counter++;
					isinList = false;
					for (auto ch_loc : bn->users_vec[i]->checkin_locations) {
						if (poi_index == ch_loc.first) {
							isinList = true;
							break;
						}

					}
					//if it seached for 1000 times and didn't find any neighbors then pass
					if (counter > 1000) {
						availble_euclidean_dist_btwn_user_POI = false;
						break;
					}
				}
				if (!availble_euclidean_dist_btwn_user_POI)
					break;
			}
			else {
				first_poi_index = poi_index;
			}
			int frequncy = 0;
			while (frequncy <= 0)
				frequncy = isUniformFrequency ? random1::uniform_int_pos(1, max_frequncy) : random1::gaussian_pos(1, max_frequncy, max_frequncy / 2.0, max_frequncy / 5.0);// assign random frequency to checkin_loc
			bn->users_vec[i]->checkin_locations[poi_index] = frequncy;
		}
	}
}
void SampleData::fix_dblp_ungraph()
{
	for (int i = 0; i < bn->no_of_users; i++) {
		bn->users_vec.push_back(new User(i, 5, toatl_n_of_used_keys));
	}
	std::ifstream fileNode("dataset/sn/com_dblp_ungraph.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file!\n";
		return;
	}
	std::string delimiter = " ";
	std::string line;
	std::map<int, int> nodes;//id, index
	while (getline(fileNode, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));
		nodes[id1] = 0;
		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		nodes[id2] = 0;
	}
	fileNode.close();

	int i = 0;
	for (auto n : nodes) {
		nodes[n.first] = i++;
	}
	std::string nodes_to_print = "";

	std::ifstream fileNode2("dataset/sn/com_dblp_ungraph.txt");

	if (!fileNode2.is_open())
	{
		std::cout << "Could not open file!\n";
		return;
	}
	while (getline(fileNode2, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		nodes_to_print += std::to_string(nodes[id1]) + " " + std::to_string(nodes[id2]) + "\n";
	}
	fileNode2.close();

	Print_To_File("com_dblp_ungraph_by_index.txt", nodes_to_print);
}
double SampleData::Load_dblp_ungraph_add_info(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequncy, bool isUniformFrequency)
{
	Generate_Social_Network_Random_Users(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequncy, isUniformFrequency);

	std::ifstream fileNode("dataset/sn/com_dblp_ungraph_by_index.txt");

	if (!fileNode.is_open())
	{
		std::cout << "Could not open file!\n";
		return 0;
	}
	std::string delimiter = " ";
	std::string line;
	while (getline(fileNode, line))
	{
		int id1 = std::stoi(line.substr(0, line.find(delimiter)));
		line.erase(0, line.find(delimiter) + 1);
		int id2 = std::stoi(line.substr(0, line.find(delimiter)));
		float random_influnce = random1::gaussian_pos(1, 10, 10 / 2.0, 10) / 10;
		bn->users_vec[id1]->outNeighbors[id2] = random_influnce;
		bn->users_vec[id2]->inNeighbors[id1] = random_influnce;
		bn->Set_edge_sup(id1, id2, 0);
	}
	fileNode.close();
	//Print_SN(true, "email-Eu-core_sn");

	return 0;
}
bool has_edge(int u, int v, const std::map<std::pair<int, int>, float>& checked_edges) {
	auto edge = std::make_pair(u, v);
	if (checked_edges.size() == 0)
		return false;
	return checked_edges.count(edge);
}
void SampleData::generate_Socail_Netword_with_k_d_truss(int max_no_checkin_loc, int max_euclidean_dist_btwn_user_POI, int max_frequency, bool isUniformFrequency, std::string filename, std::vector<int> users_size_vec, int _min_no_follower, int _max_no_follower, int distribution_type)  // uniform=1, guassian =2 , zipf=3
{
	std::string file_name = std::to_string(total_n_of_intersection_points) + "_nodes.txt";
	std::ifstream fileRnodes("dataset/synthatic/rn/" + file_name);
	if (!fileRnodes.is_open())
	{
		std::cout << "Could not open file!\n";
	}
	std::string line;
	int i;
	while (getline(fileRnodes, line))
	{
		i = std::stoi(line.substr(0, line.find(" ")));
		line.erase(0, line.find(" ") + 1);
		intersection_points.at(2 * i) = std::stof(line.substr(0, line.find(" ")));
		line.erase(0, line.find(" ") + 1);
		intersection_points.at(2 * i + 1) = std::stof(line.substr(0, line.find(" ")));
	}
	fileRnodes.close();


	file_name = std::to_string(total_n_of_intersection_points) + "_edges.txt";
	std::ifstream fileRedgest("dataset/synthatic/rn/" + file_name);
	if (!fileRedgest.is_open())
	{
		std::cout << "Could not open file!\n";
	}

	int st_node, end_node;
	float dist;
	while (getline(fileRedgest, line))
	{
		st_node = std::stoi(line.substr(0, line.find("-")));
		line.erase(0, line.find("-") + 1);
		end_node = std::stoi(line.substr(0, line.find("=")));
		line.erase(0, line.find("=") + 1);
		dist = std::stof(line.substr(0, line.find(" ")));
		roadNetworkAdjMtrx[st_node].insert(std::make_pair(end_node, dist));
	}
	fileRedgest.close();

	std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";
	Load_POI("dataset/poi/" + distribution_type_name + "ip_" + std::to_string(total_n_of_intersection_points) + "_poi_" + std::to_string(no_of_POI) + ".txt");



	uint64_t seed = 100;
	int k = 6;
	SocialNetwork net(seed, k, _min_no_follower, _max_no_follower);
	// random community diameters in [2,6]
	net.init_communities(20, 2, 6);

	std::map<std::pair<int, int>, float> checked_edges;
	for (auto users_size : users_size_vec) {
		int start_index = bn->users_vec.size();
		net.expand_to(users_size);

		Generate_Social_Network_Random_Users_add_extra(users_size);
		for (int u = 0; u < net.G.n; u++)
			for (int v : net.G.adj[u]) {
				if (u >= v) continue;
				auto edge = std::make_pair(u, v);
				if (has_edge(u, v, checked_edges)) {
					float random_influnce = checked_edges[edge];
					bn->users_vec[u]->outNeighbors[v] = random_influnce;
					bn->users_vec[v]->inNeighbors[u] = random_influnce;
					bn->Set_edge_sup(u, v, 0);
				}
				else {
					float random_influnce = random1::gaussian_pos(1, 10, 10 / 2.0, 10) / 10;
					bn->users_vec[u]->outNeighbors[v] = random_influnce;
					bn->users_vec[v]->inNeighbors[u] = random_influnce;
					bn->Set_edge_sup(u, v, 0);
					checked_edges[edge] = random_influnce;
				}
			}
		Assign_poi_to_users(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, isUniformFrequency, users_size, start_index);
		Print_SN(true, filename);
	}
}

void SampleData::Print_SN_timestamp(bool IsWrite, std::string filename)
{
	printf("\SN\n");
	if (!IsWrite) {
		filename = filename + "_sn.txt";
		std::string text = bn->print();
		std::ofstream myfile;
		myfile.open(filename);
		myfile << text;
		myfile.close();
	}
	else {
		std::string file_name = filename + "ip_" + std::to_string(total_n_of_intersection_points) + "_users_" + std::to_string(bn->no_of_users) + "_poi_" + std::to_string(no_of_POI) + "_sn_timestamp.txt";
		std::ofstream myfile;
		myfile.open(file_name);
		for (int i = 0; i < bn->no_of_users; i++) {
			User* u = bn->users_vec[i];
			myfile << u->id;
			for (auto& u_ch_in : u->checkin_locations) {
				myfile << " " << u_ch_in.first << ":" << u_ch_in.second;
			}
			myfile << "|";
			for (auto& u_out_n : u->timed_checkins) {
				myfile << " " << u_out_n.first << ":";
				for (auto& t : u_out_n.second)
					myfile << t << ",";
			}
			myfile << "\n";
		}
		myfile.close();
	}
}
