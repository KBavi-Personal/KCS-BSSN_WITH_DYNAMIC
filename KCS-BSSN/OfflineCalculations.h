#pragma once
#include "SampleData.h"
#include "rand.h"
class OfflineCalculations
{
public:
	SampleData* sd;


	OfflineCalculations(SampleData *_sd);
	~OfflineCalculations();
	void Run(int pivot_iterations, std::string filename, int distribution_type);
	int Get_ub_f_sum(User* u);
	float Get_ub_f_avg(User* u);
	float Get_ub_w_in(User* u);
	float Get_ub_w_out(User* u);
	float Get_ub_sub(User* u);
	float Get_largest_f_avg();

	void SN_Pivot_Selection(int iteration);

	void RN_Pivot_Selection(int iteration);

	void Calculate_edges_sup_in_SN2();

	void Print_SN_Pivot(std::string filename);
	void Print_RN_Pivot(std::string filename);
	void Compute_Key_fsum_fmax_for_users();
	void Compute_poi_max_freq();
	void Compute_POI_fsum_count_per_key();
private:
	std::vector<int> poi_max_freq;
};

