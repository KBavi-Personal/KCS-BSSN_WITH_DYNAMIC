#pragma once
#include <vector>
#include <set>
#include <string>
#include<unordered_set>

#include "Sextuple.h"
#include "CommStru.h"

struct SampleData;
class Dijkstra;

class Ref2 {
public:
	explicit Ref2(SampleData* _sd);
	~Ref2();

	std::vector<int> dist_tmp;
	// core
	void Run_Ref2(const std::vector<int>& Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index, Sextuple<bool, bool, bool, bool, bool, bool, bool> m);

	// Ref2s
	void Keyword_Based_Ref2(const std::vector<int>& Q);
	void Omega_Based_Ref2(int omega);
	void Pi_Based_Ref2(float pi);
	void Social_Distance_Based_Ref2(int d, int q_index);
	void Structural_Cohesiveness_Ref2(int k);
	void Spatial_Distance_Ref2_cache(float sigma, const std::set<int>& v_p_prime_set, int q_index);

	// influence
	bool HasInflunceToAllRemaining(Dijkstra* dijkstra, float theta, int source_index, int q_index);

	void Influence_Based_Ref2_dijkstra(Dijkstra* dijkstra, float theta, int q_index);


	void remove_user_from_remaining_users(int index);
	void remove_pruned_user();
	void remove_pruned_POI();	
	
	std::set<int> Reset();
	std::string  print_remain_candidates(std::string label);
	

	void Run_Ref_dynamic(const std::vector<int>& Q, int k, int d, int omega, float pi, float theta, float sigma, int q_index);
	void Keyword_Based_Ref2_dynamic(const std::vector<int>& Q);




private:
	SampleData* sd = nullptr;
};

