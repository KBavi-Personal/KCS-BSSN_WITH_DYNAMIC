#include "Dijkstra.h"
#include "BipartiteNetwork.h"
#include "SampleData.h"
#include "Triple.h"
#include "MinHeap.h"
#include <queue>
#include <cfloat>

Dijkstra::Dijkstra(SampleData* _sd)
{
	sd = _sd;
}

Dijkstra::~Dijkstra()
{
}
float Dijkstra::Calculate_shotest_dist_road(SampleData* sd, int src_loc_index, int destination_loc_index)
{
	MinHeap* heap = new MinHeap(sd->total_n_of_intersection_points);

	std::vector <Triple<float, int, bool>> dist_counter_isdeleted_vec;// total distance and counter
	Triple<int, float, int> heap_element = Triple<int, float, int>(0, 0.0, 0);
	for (int ip = 0; ip < sd->total_n_of_intersection_points; ip++)
	{
		if (src_loc_index == ip)
			dist_counter_isdeleted_vec.push_back(Triple<float, int, bool>(0.0, 0, false));
		else
			dist_counter_isdeleted_vec.push_back(Triple<float, int, bool>(FLT_MAX, 0, false));
		heap_element.first = ip;
		heap_element.second = dist_counter_isdeleted_vec[ip].first;
		heap->insert_minheap(heap_element);
	}
	bool IsBreak = false;
	float dist = 0.0;
	double weight = 0.0;
	while (heap->size > 0)
	{
		Triple<int, float, int> heap_top_element = heap->get_min();
		// remove the top heap------------------------------------ -
		heap->delete_minimum();
		dist_counter_isdeleted_vec[heap_top_element.first].third = true;

		std::map<int, float> neighbors = sd->roadNetworkAdjMtrx[heap_top_element.first];
		for (auto neighbor : neighbors) {
			weight = neighbor.second;
			if (!dist_counter_isdeleted_vec[neighbor.first].third) {
				int element_index = heap->find_element_index_value(neighbor.first);
				if (element_index != -1 && dist_counter_isdeleted_vec[heap_top_element.first].first != FLT_MAX && dist_counter_isdeleted_vec[heap_top_element.first].first + weight < dist_counter_isdeleted_vec[neighbor.first].first)
				{
					dist_counter_isdeleted_vec[neighbor.first].first = dist_counter_isdeleted_vec[heap_top_element.first].first + weight;
					heap->arr[element_index].second = dist_counter_isdeleted_vec[neighbor.first].first;
					heap->update_element(element_index);
				}
			}
			dist_counter_isdeleted_vec[neighbor.first].second = dist_counter_isdeleted_vec[neighbor.first].second + 1;
			if (destination_loc_index == neighbor.first && sd->roadNetworkAdjMtrx[neighbor.first].size() == dist_counter_isdeleted_vec[neighbor.first].second) {
				IsBreak = true;
				dist = dist_counter_isdeleted_vec[neighbor.first].first;
				break;
			}
		}
		if (IsBreak) break;
	}
	delete heap;
	return dist;
}
void Dijkstra::Calculate_shotest_dist_road_s_to_all(int pr_index)
{
	if (pr_index < 0 || pr_index >= (int)sd->P_r_vec.size()) return;

	const int N = sd->total_n_of_intersection_points;
	const int src = sd->P_r_vec[pr_index];
	if (src < 0 || src >= N) return;

	const float INF = std::numeric_limits<float>::infinity();

	// dist to all intersections
	std::vector<float> dist(N, INF);
	dist[src] = 0.0f;

	using Node = std::pair<float, int>; // (dist, node)
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
	pq.push({ 0.0f, src });

	while (!pq.empty()) {
		auto [du, u] = pq.top();
		pq.pop();

		if (du != dist[u]) continue; // stale entry

		const auto& neighbors = sd->roadNetworkAdjMtrx[u];

		for (const auto& [v, w] : neighbors) {
			const float nd = du + w;
			if (nd < dist[v]) {
				dist[v] = nd;
				pq.push({ nd, v });
			}
		}
	}


	for (int i = 0; i < sd->no_of_POI; ++i) {
		const int ip = sd->poi_vec[i]->my_closest_intesect_point;
		sd->poi_vec[i]->dist_P_r[pr_index] = (ip >= 0 && ip < N) ? dist[ip] : INF;
	}

}
void Dijkstra::Calculate_shotest_dist_road_All_intersection_points(SampleData* sd)
{
	std::vector <Triple<float, int, bool>> dist_counter_isdeleted_vec;// total distance and counter
	Triple<int, float, int> heap_element = Triple<int, float, int>(0, 0.0, 0);
	for (int intersection_point = 0; intersection_point < 1000; intersection_point++) {

		dist_counter_isdeleted_vec.clear();
		MinHeap* heap = new MinHeap(sd->total_n_of_intersection_points);



		for (int ip = 0; ip < sd->total_n_of_intersection_points; ip++)
		{
			if (intersection_point == ip)
				dist_counter_isdeleted_vec.push_back(Triple<float, int, bool>(0.0, 0, false));
			else
				dist_counter_isdeleted_vec.push_back(Triple<float, int, bool>(FLT_MAX, 0, false));
			if (ip < intersection_point)continue;
			heap_element.first = ip;
			heap_element.second = dist_counter_isdeleted_vec[ip].first;
			heap->insert_minheap(heap_element);
		}
		bool IsBreak = false;
		float dist = 0.0;
		double weight = 0.0;
		while (heap->size > 0)
		{
			Triple<int, float, int> heap_top_element = heap->get_min();

			heap->delete_minimum();
			dist_counter_isdeleted_vec[heap_top_element.first].third = true;

			std::map<int, float> neighbors = sd->roadNetworkAdjMtrx[heap_top_element.first];
			for (std::pair<int, float> neighbor : neighbors) {
				weight = neighbor.second;

				if (!dist_counter_isdeleted_vec[neighbor.first].third) {
					int element_index = heap->find_element_index_value(neighbor.first);
					if (element_index != -1 && dist_counter_isdeleted_vec[heap_top_element.first].first != FLT_MAX && dist_counter_isdeleted_vec[heap_top_element.first].first + weight < dist_counter_isdeleted_vec[neighbor.first].first)
					{
						dist_counter_isdeleted_vec[neighbor.first].first = dist_counter_isdeleted_vec[heap_top_element.first].first + weight;
						heap->arr[element_index].second = dist_counter_isdeleted_vec[neighbor.first].first;
						heap->update_element(element_index);
					}
				}
				dist_counter_isdeleted_vec[neighbor.first].second = dist_counter_isdeleted_vec[neighbor.first].second + 1;
			}
		}
		delete heap;
	}
}

std::vector<float> Dijkstra::MaxInfluenceFromSource_All_OutNeighbors(int src_user_Id, float theta)
{
	const int N = sd->bn->no_of_users;

	std::vector<int> idx(N, -1);
	idx.reserve(N);
	for (int i = 0; i < sd->no_remining_users; ++i) {
		idx[sd->remining_users[i]] = i;
	}

	std::vector<float> best(N, 0.0f); // influence, 0 if unreachable / pruned
	if (src_user_Id < 0 || src_user_Id >= N) return best;
	if (idx[src_user_Id] == -1) return best;

	const float INF = std::numeric_limits<float>::infinity();
	std::vector<float> dist(sd->no_remining_users, INF);
	dist[idx[src_user_Id]] = -1.0f;

	using Node = std::pair<float, int>; // (dist, userId)
	std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
	pq.push({ -1.0f, src_user_Id });

	while (!pq.empty()) {
		auto [du, uId] = pq.top();
		pq.pop();

		const int uIdx = idx[uId];
		if (uIdx == -1) continue;
		if (du != dist[uIdx]) continue;

		const float inf_u = -du;
		if (inf_u < theta) break;

		best[uId] = inf_u;

		for (const auto& nb : sd->bn->users_vec[uId]->outNeighbors) {
			const int vId = nb.first;
			const float w = (float)nb.second;

			if (vId < 0 || vId >= N) continue;
			if (sd->bn->users_vec[vId]->isPruned) continue;
			const int vIdx = idx[vId];
			if (vIdx == -1) continue;

			const float nd = du * w; // du is negative
			if (nd < dist[vIdx]) {
				dist[vIdx] = nd;
				pq.push({ nd, vId });
			}
		}
	}

	return best;
}



