#include <iostream>
#include <vector>
#include "SampleData.h"
#include "BipartiteNetwork.h"
#include "Dijkstra.h"
#include "Indexing.h"
#include "OfflineCalculations.h"
#include "BFS.h"
#include <chrono>
#include "ref2.h"
#include "SocialNetworkGraphKTruss.h"
#include <sstream>
#include "Tree.h"
#include <queue>
#include "Pruning.h"
#include "TreeTextIO.h"
#include "CommStru.h"
#include "DynamicApproch.h"

long double index_time[10], ind_st, ind_ed;


#define BASE_LINE_TIME 0
#define KCSBSSN_TIME 1
#define INDEXING_TIME 2
#define OFFLINE_TIME 3
#define FIND_SHORTEST_R_DEST_TIME 4
#define FIND_SHORTEST_S_DEST_TIME 5
#define GENERATE_BN_TIME 6
std::bitset<500000> GLOBAL_CHECKED_q;
std::vector<int> GLOBAL_TESTED_POS_q;
int GLOBAL_CHECHED_TESTED_q_COUNTER;

// Parameter values
std::vector<int> k{ 2, 3, 4, 5,6 };
std::vector<int> d{ 1, 2, 3, 4, 5 };
std::vector<float> omega{ 0.2f,0.4f, 0.6f, 0.7f, 0.9f };
std::vector<float> pi{ 0.2f,0.4f, 0.6f, 0.7f, 0.9f };
std::vector<float> theta{ 0.2f,0.4f, 0.6f, 0.7f, 0.9f };
//std::vector<float> omega{ 0.1f, 0.2f,0.3f,0.4f, 0.5f,0.6f, 0.7f,0.8f, 0.9f,1 };
//std::vector<float> pi{ 0.1f, 0.2f,0.3f,0.4f, 0.5f,0.6f, 0.7f,0.8f, 0.9f,1 };
//std::vector<float> theta{ 0.1f, 0.2f,0.3f,0.4f, 0.5f,0.6f, 0.7f,0.8f, 0.9f,1 };
std::vector<float> sigma{ 1,2,3, 4,5,6 };




int Qset_size = 5;
int k_size = k.size();
int d_size = d.size();
int omega_size = omega.size();
int pi_size = pi.size();
int theta_size = theta.size();
int sigma_size = sigma.size();



int default_Q = 2;
int default_k = 1;
int default_d = 2;
int default_omega = 3;
int default_pi = 3;
int default_theta = 3;
int default_sigma = 5;


bool IsItemInVec(std::vector<int> vec, int item) {
	for (auto i : vec)
		if (i == item)return true;
	return false;
}
std::vector<std::vector<int>>  Get_Q_sets(int Q_size, int total_no_of_keys) {
	std::vector<int> Q_set_main;
	std::vector<std::vector<int>> Q_sets;
	bool cal_rand_Q;
	std::cout << "Get random Q set (1)-yes, (0)-no" << std::endl;
	std::cin >> cal_rand_Q;

	if (cal_rand_Q) {
		for (int i = 0; i < Q_size; i++) {
			int key = random1::uniform_int_pos(0, total_no_of_keys);
			while (IsItemInVec(Q_set_main, key))
				key = random1::uniform_int_pos(0, total_no_of_keys);
			Q_set_main.push_back(key);
		}
	}
	else {
		int key;
		std::cout << "Enter " << Q_size << " keys" << std::endl;
		for (int i = 0; i < Q_size; i++) {
			std::cin >> key;
			Q_set_main.push_back(key);
		}
	}
	if (cal_rand_Q) {
		std::cout << "Q_set= ";
		for (int i = 0; i < Q_size; i++) {
			std::cout << Q_set_main[i] << " ";
		}
		std::cout << "\n";
	}
	std::vector<int> Q_set;
	Q_set.push_back(Q_set_main[0]);
	Q_set.push_back(Q_set_main[1]);
	Q_set.push_back(Q_set_main[2]);
	Q_sets.push_back(Q_set);
	Q_set.push_back(Q_set_main[3]);
	Q_set.push_back(Q_set_main[4]);
	Q_sets.push_back(Q_set);
	Q_set.push_back(Q_set_main[5]);
	Q_set.push_back(Q_set_main[6]);
	Q_sets.push_back(Q_set);
	Q_set.push_back(Q_set_main[7]);
	Q_set.push_back(Q_set_main[8]);
	Q_sets.push_back(Q_set);
	Q_set.push_back(Q_set_main[9]);
	Q_set.push_back(Q_set_main[10]);
	Q_sets.push_back(Q_set);
	/*std::vector<int> Q_set1;
	for (int i = 0; i < 50; i++) {

		Q_set1.push_back(i);

	}
	Q_sets.push_back(Q_set1);*/

	/*for (int i = 1; i <= Q_size / 2; i++) {
		std::vector<int> Q_set;
		for (int j = 0; j < i * 2; j += 2) {
			Q_set.push_back(Q_set_main[j]);
			Q_set.push_back(Q_set_main[j + 1]);
		}
		Q_sets.push_back(Q_set);
	}*/
	return Q_sets;
}
std::set<int> Get_q_random_community(SampleData* sd, int q_index, int d)
{
	std::vector<int> dist_tmp;

	std::set<int> community;
	float random_num = 0;
	// Create a queue for BFS
	std::queue<std::pair<int, int>> queue;//user_id, level

	std::pair<int, int> s;
	queue.push(std::make_pair(q_index, 0));
	while (!queue.empty()) {
		// Dequeue a vertex from queue 
		s = queue.front();
		if (s.second >= d)break;
		queue.pop();

		BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, s.first, 1, dist_tmp, false);
		for (int i = 0; i < sd->bn->no_of_users; ++i) {
			if (dist_tmp[i] == -1 || i == s.first) continue;
			random_num = random1::uniform_int_pos(0, 10);
			if (random_num >= 5) {
				if (community.count(i) == 1) continue;
				community.insert(i);
				queue.push(std::make_pair(i, s.second + 1));
			}
		}
	}
	return community;
}
void GenerateOfflineData(SampleData* sd, int pivot_iterations, std::string filename, int distribution_type) {

	//sd->Print_To_File(filename, "#users=" + std::to_string(sd->bn->no_of_users) + ", #POI= " + std::to_string(sd->no_of_POI) + ", #Intersections= " + std::to_string(sd->total_n_of_intersection_points));
	double spent_time = 0.0;

	OfflineCalculations* offlineCal = new OfflineCalculations(sd);
	offlineCal->Run(pivot_iterations, filename, distribution_type);
	delete offlineCal;
}
int Select_q_index(int n_users) {
	int index;
	if (GLOBAL_CHECHED_TESTED_q_COUNTER < GLOBAL_TESTED_POS_q.size()) {
		index = GLOBAL_TESTED_POS_q[GLOBAL_CHECHED_TESTED_q_COUNTER];
		GLOBAL_CHECHED_TESTED_q_COUNTER++;
	}
	else {
		index = random1::uniform_int_pos(0, n_users - 1);
	}
	while (GLOBAL_CHECKED_q.test(index)) {
		index = random1::uniform_int_pos(0, n_users - 1);
	}
	GLOBAL_CHECKED_q[index] = true;
	return index;
}
std::string RunPrametersPowerEva(SampleData* sd, Tree* tree, std::vector<int> Q, int k, int d, float omega, float pi, float theta, float sigma, std::string filename, Ref2* ref, int distribution_type, Sextuple<bool, bool, bool, bool, bool, bool, bool> m) {
	int omega_int = omega * sd->larget_f_sum;
	float pi_fl = pi * sd->larget_f_avg;
	if (pi_fl < 1)pi_fl = 1;


	std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";

	int total_remain_users = 0;
	int Run_ref2_user = 0;
	std::vector<int> dist_tmp;

	int q_tested = 0;
	int text_all = 0;
	int q_index;

	std::ostringstream results_to_print;
	int counter = 0;

	q_index = -1;
	while (q_index < sd->bn->no_of_users - 1 && q_tested < 5) {

		q_index++;
		counter++;

		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		sd->Reset_POI_pruned_flag();
		sd->bn->Reset_Users_pruned_flag();
		for (auto u_id : sd->bn->users_vec) {
			sd->remining_users[sd->no_remining_users] = u_id->id;
			sd->no_remining_users++;
		}
		int remain_users2 = 0;
		int remian_poi = 0;

		tree->KCS_BSSN_Query_Answer(Q, k, d, omega_int, pi_fl, theta, sigma, q_index, m);
		ref->Run_Ref2(Q, k, d, omega_int, pi_fl, theta, sigma, q_index, m);

		remain_users2 = sd->no_remining_users;
		remian_poi = sd->no_remining_POI;
		if (sd->no_remining_users < 3) continue;
		q_tested++;

		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		sd->Reset_POI_pruned_flag();
		sd->bn->Reset_Users_pruned_flag();
		for (auto u_id : sd->bn->users_vec) {
			sd->remining_users[sd->no_remining_users] = u_id->id;
			sd->no_remining_users++;
		}
		results_to_print << distribution_type_name << filename << "," << sd->bn->no_of_users << "," << Q.size() << "," << k << "," << d << "," << omega_int << "," << omega << "," << pi_fl << "," << pi << "," << theta << "," << sigma << "," << q_index << "," << "," << remain_users2 << "," << remian_poi << "\n";
		total_remain_users += remain_users2;


	}
	if (q_tested == 0) q_tested = 1;
	std::ostringstream  m_options;
	m_options << (m.keyword ? "keyword " : "") << (m.pi ? "pi " : "") << (m.omega ? "omega " : "") << (m.social ? "social " : "") << (m.structrual ? "structural " : "") << (m.influnce ? "influence " : "") << (m.spacial ? "spatial " : "");
	double avg_remain_users = total_remain_users / q_tested;
	results_to_print << distribution_type_name << filename << "," << sd->bn->no_of_users << "," << Q.size() << "," << k << "," << d << "," << omega_int << "," << omega << "," << pi_fl << "," << pi << "," << theta << "," << sigma <<
		",avg," << "," << avg_remain_users << "," << m_options.str() << "\n";
	return results_to_print.str();

}
std::string Run(SampleData* sd, Tree* tree, std::vector<int> Q, int k, int d, float omega, float pi, float theta, float sigma, std::string filename, Ref2* ref, bool isDefualts, int distribution_type) {
	int omega_int = omega * sd->larget_f_sum;
	float pi_fl = pi * sd->larget_f_avg;
	if (pi_fl < 1)pi_fl = 1;

	std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";

	GLOBAL_CHECKED_q.reset();
	GLOBAL_CHECHED_TESTED_q_COUNTER = 0;


	double basline_time_total = 0.0;
	double KCS_BSSN_Query_time_total = 0.0;
	double refinement_time_total = 0.0;
	double basline_time_for_all_comm = 0.0;
	double total_basline_time_for_all_comm = 0.0;
	int total_users_up_to_d_distance = 0;

	int users_up_to_d_distance = 0;
	int Run_ref2_user = 0; double Run_ref2_time = 0;
	int tested_q_for_basline = 0;
	std::vector<int> dist_tmp;

	int q_tested = 0;
	int text_all = 0;
	int q_index;

	std::ostringstream results_to_print;
	std::ostringstream temp;
	int counter = 0;
	int max_try = 9000;
	Sextuple m = { true, true, true, true,true,true,true };

	while (q_tested < 5) {

		counter++;
		q_index = Select_q_index(sd->bn->no_of_users - 1);

		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		sd->Reset_POI_pruned_flag();
		sd->bn->Reset_Users_pruned_flag();
		for (auto u_id : sd->bn->users_vec) {
			sd->remining_users[sd->no_remining_users] = u_id->id;
			sd->no_remining_users++;
		}
		double KCS_BSSN_Query_time = 0.0;
		const auto start5 = std::chrono::high_resolution_clock::now();
		tree->KCS_BSSN_Query_Answer(Q, k, d, omega_int, pi_fl, theta, sigma, q_index, m);
		ref->Run_Ref2(Q, k, d, omega_int, pi_fl, theta, sigma, q_index, m);
		const auto end5 = std::chrono::high_resolution_clock::now();
		const std::chrono::duration<double, std::milli> duration_ms5 = end5 - start5;
		KCS_BSSN_Query_time = duration_ms5.count();


		//pass after max_try attempts
		if (sd->no_remining_users == 0 && counter < max_try) continue;
		if (sd->no_remining_users > 0) GLOBAL_TESTED_POS_q.push_back(q_index);

		q_tested++;
		bool is_in_community = false;
		double basline_time = 0.0;
		int num_tested_samples = 5;
		int num_tested_samples_not_zero = 0;
		if (isDefualts) {
			for (int j = 0; j < num_tested_samples; j++) {
				const auto start1 = std::chrono::high_resolution_clock::now();
				std::set<int> community = Get_q_random_community(sd, q_index, d);
				sd->no_remining_users = 0;
				sd->no_remining_POI = 0;
				sd->Reset_POI_pruned_flag();
				sd->bn->Reset_Users_pruned_flag();
				if (community.size() > 2)
				{
					num_tested_samples_not_zero++;
					for (int k = 0; k < sd->bn->no_of_users; k++)
						sd->bn->users_vec[k]->isPruned = community.count(k) == 0;
					for (auto u_id : community) {
						sd->remining_users[sd->no_remining_users] = u_id;
						sd->no_remining_users++;
					}

					ref->Run_Ref2(Q, k, d, omega_int, pi_fl, theta, sigma, q_index, m);
					const auto end1 = std::chrono::high_resolution_clock::now();
					const std::chrono::duration<double, std::milli> duration_ms1 = end1 - start1;
					basline_time += duration_ms1.count();
				}
			}
		}



		sd->no_remining_users = 0;
		sd->no_remining_POI = 0;
		sd->Reset_POI_pruned_flag();
		sd->bn->Reset_Users_pruned_flag();
		for (auto u_id : sd->bn->users_vec) {
			sd->remining_users[sd->no_remining_users] = u_id->id;
			sd->no_remining_users++;
		}

		users_up_to_d_distance = 0;
		BFS::ComputeDistancesFromSourceCutoff(sd->bn->users_vec, q_index, d, dist_tmp, false);
		for (int i = 0; i < sd->bn->no_of_users; ++i) {
			if (dist_tmp[i] != 0) users_up_to_d_distance++;
		}

		total_users_up_to_d_distance += users_up_to_d_distance;

		if (num_tested_samples_not_zero > 0) {
			basline_time = basline_time / num_tested_samples_not_zero;
			tested_q_for_basline++;
		}

		results_to_print << distribution_type_name << filename << "," << sd->bn->no_of_users << "," << Q.size() << "," << k << "," << d << "," << omega_int << "," << omega << "," << pi_fl << "," << pi << "," << theta << "," << sigma << "," << q_index << "," << basline_time << "," << users_up_to_d_distance << "," << KCS_BSSN_Query_time << "\n";

		basline_time_total += basline_time;
		KCS_BSSN_Query_time_total += KCS_BSSN_Query_time;
		total_basline_time_for_all_comm += basline_time_for_all_comm;


	}
	if (q_tested == 0) q_tested = 1;


	double avg_basline_time = tested_q_for_basline > 0 ? basline_time_total / tested_q_for_basline : 0;
	double avg_KCS_BSSN_time = KCS_BSSN_Query_time_total / q_tested;
	double avg_users_up_to_d_distance = total_users_up_to_d_distance / q_tested;

	//"type,user_n0,,Q,k,d,omega,omega_nor,pi,pi_nor,theta,sigma,q,Basline_time_one_comm,users_up_to_d_distance,KCS_BSSN_Query_time,";
	results_to_print << distribution_type_name << filename << "," << sd->bn->no_of_users << "," << Q.size() << "," << k << "," << d << "," << omega_int << "," << omega << "," << pi_fl << "," << pi << "," << theta << "," << sigma <<
		",avg," << avg_basline_time << "," << avg_users_up_to_d_distance << "," << avg_KCS_BSSN_time << "\n";

	return results_to_print.str();

}
Tree* GenerateTree(SampleData* sd, int pivot_iterations, int pivot_size, std::string filename, int level_dev, float bs_w, float ss_w, float rs_w) {

	Tree* tree = new Tree(sd);
	tree->BuildTree(pivot_iterations, level_dev, pivot_size, bs_w, ss_w, rs_w);
	return tree;
}
void RunAllPossibleParameteres(SampleData* sd, std::vector<std::vector<int>>& Q_sets, std::string filename, bool justdefault, bool is_pruning_pow_eva, int distribution_type, bool loadTree) {


	std::cout << filename + "_" << sd->bn->no_of_users << "_" << sd->no_of_POI << std::endl;
	int pivot_iterations = 5;
	int tree_piv_size[] = { 100 };
	Tree* tree = NULL;
	Ref2* ref = new Ref2(sd);

	double offline_time = 0.0;
	const auto start = std::chrono::high_resolution_clock::now();
	GenerateOfflineData(sd, pivot_iterations, filename, distribution_type);
	const auto end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms = end - start;
	offline_time = duration_ms.count();

	std::string tree_filename;
	float bs_w = 1, ss_w = 2, rs_w = 1;
	int tree_piv_size_single = sd->bn->no_of_users / 100, level_dev = 2;
	pivot_iterations = 5;

	std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";
	tree_filename = distribution_type_name + filename + "_" + std::to_string(tree_piv_size_single) + "piv_result";
	loadTree = false;

	double GenerateTree_time = 0.0;
	const auto start1 = std::chrono::high_resolution_clock::now();
	tree = GenerateTree(sd, pivot_iterations, tree_piv_size_single, filename, level_dev, bs_w, ss_w, rs_w);
	const auto end1 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> duration_ms1 = end1 - start1;
	GenerateTree_time = duration_ms1.count();

	std::ostringstream  time_lbl;
	time_lbl << distribution_type_name << sd->bn->no_of_users << "," << offline_time / 1000 << "," << GenerateTree_time / 1000 << "," << (GenerateTree_time + offline_time) / 1000 << "\n";
	sd->Print_To_File(filename + "_offline_time.txt", time_lbl.str());


	std::cout << "Query is running..." << std::endl;

	std::string results_to_print_titles = "type,user_no,Q,k,d,omega,omega_nor,pi,pi_nor,theta,sigma,q,Basline_time_one_comm,users_up_to_d_distance,KCS_BSSN_Query_time";

	std::string results_to_print = results_to_print_titles;
	results_to_print += " \n";


	Sextuple m = { true, true, true, true,true,true,true };
	std::string  ll;
	bool change_defaults = false;

	if (is_pruning_pow_eva) {

		m = { true, false, false, false,false,false,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, false, false,false,false,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, true, false,false,false,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, true, true,false,false,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, true, true,true,false,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, true, true,true,true,false };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);
		m = { true, true, true, true,true,true,true };
		results_to_print += RunPrametersPowerEva(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, distribution_type, m);

		sd->Print_To_File(filename + ll + "_prouningPower.csv", results_to_print);
		results_to_print = "";
		ll = "";
	}
	else {
		if (justdefault) {
			results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, true, distribution_type);
		}
		else {
			results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, true, distribution_type);
			for (int i1 = 0; i1 < Qset_size; i1++)
				results_to_print += Run(sd, tree, Q_sets[i1], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i2 = 0; i2 < k_size; i2++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[i2], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i3 = 0; i3 < d_size; i3++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[i3], omega[default_omega], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i4 = 0; i4 < omega_size; i4++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[i4], pi[default_pi], theta[default_theta], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i5 = 0; i5 < pi_size; i5++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[i5], theta[default_theta], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i6 = 0; i6 < theta_size; i6++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[i6], sigma[default_sigma], filename, ref, false, distribution_type);
			for (int i7 = 0; i7 < sigma_size; i7++)
				results_to_print += Run(sd, tree, Q_sets[default_Q], k[default_k], d[default_d], omega[default_omega], pi[default_pi], theta[default_theta], sigma[i7], filename, ref, false, distribution_type);
		}
		sd->Print_To_File(filename + ll + "_result.csv", results_to_print);
		results_to_print = "";
		ll = "";
	}


	delete tree;
	delete ref;
}
std::map<std::pair<int, int>, float> checked_edges;
void Change_Influnce(SampleData* sd, std::string filename)
{
	float random_influnce = 0;
	for (int u = 0; u < sd->bn->no_of_users; u++) {
		for (auto e : sd->bn->users_vec[u]->outNeighbors) {
			int v = e.first;

			auto edge = std::make_pair(u, v);

			if (checked_edges.size() > 0 && checked_edges.count({ u, v }))
				random_influnce = checked_edges[{u, v}];
			else
				random_influnce = random1::gaussian_pos(5, 10, 5.0, 10) / 10;
			sd->bn->users_vec[u]->outNeighbors[v] = random_influnce;
			sd->bn->users_vec[v]->inNeighbors[u] = random_influnce;
			checked_edges[{u, v}] = random_influnce;
		}
	}
	sd->Print_SN(true, filename);
}
void Print_To_File(std::string filename, std::string text)
{
	std::string file_name = filename;

	std::ofstream myfile;
	myfile.open(file_name, std::ios_base::app);
	myfile << text << "\n";
	myfile.close();

}

std::vector<int> Pick_random_users_for_dynamic(int num_users_to_pick, int  no_of_users) {
	std::vector<int> picked_users;
	//int num_users_to_pick = static_cast<int>(no_of_users * percentage);
	for (int i = 0; i < num_users_to_pick; i++) {
		int random_user_id = random1::uniform_int_pos(0, no_of_users - 1);
		picked_users.push_back(random_user_id);
	}
	return picked_users;
}
std::vector<int> Pick_random_users_for_dynamic_percentage(double percentage, int  no_of_users) {
	std::vector<int> picked_users;
	int num_users_to_pick = static_cast<int>(no_of_users * percentage);
	for (int i = 0; i < num_users_to_pick; i++) {
		int random_user_id = random1::uniform_int_pos(0, no_of_users - 1);
		picked_users.push_back(random_user_id);
	}
	return picked_users;
}

void Print_To_File_main(std::string filename, std::string text)
{
	std::string file_name = filename;

	std::ofstream myfile;
	myfile.open(file_name, std::ios_base::app);
	myfile << text << "\n";
	myfile.close();

}

int main()
{
	std::srand(std::time(0));

	SampleData* newSampleData;
	bool isReal = false;
	bool isTwitter = true;

	int no_p_r = 5;
	int no_p_s = 5;
	int _max_x = 11;
	int _max_y = 11;

	int max_no_of_keys_per_p = 8;
	int total_no_of_keys = 50;
	int _max_no_follower = 100;
	int max_no_checkin_loc = 10;
	int max_euclidean_dist_btwn_user_POI = 3;//max dist is 11
	int max_frequency = 10;


	int Q_size = 11;

	std::vector<int> number_users = { 10000,20000,30000,40000,50000,100000,200000 };

	int number_intersection_points[] = { 10000,20000,30000,40000,50000 };
	int number_POI[] = { 2000,4000,6000,8000,10000 };

	int default_number_users_index = 2;
	int default_number_Intersection_points_index = 1;
	int default_number_POI_index = 4;

	std::vector<int> Q_set_main;
	std::vector<std::vector<int>> Q_sets;
	Q_sets = Get_Q_sets(Q_size, total_no_of_keys);



	int choose = 0;
	std::cout << "Which one you want to test? Type: (1)-real_datasets,  (2)-Generate POI (3)-generate synthetic and assign POI (BiPartite) (4)-run synthetic (5)-run dynamic updates" << std::endl;
	std::cin >> choose;


	bool is_pruning_pow_eva = 0;
	std::cout << "Is this a seprate pruning power evaluation? (1)-yes, (0)-no" << std::endl;
	std::cin >> is_pruning_pow_eva;

	std::string filename;
	int real_choose = 0;
	bool load_tree = true;
	switch (choose) {
	case 1:
		total_no_of_keys = 72;
		std::cout << "Which one you want to test? Type: (1)-Twitter, (2)-Epinions (3)-dblp_un_graph" << std::endl;
		std::cin >> real_choose;
		switch (real_choose) {
		case 1:
			try {
				filename = "Twitter";
				newSampleData = new SampleData(81306, 104770, no_p_r, no_p_s, _max_x, _max_y, 21048, total_no_of_keys);
				newSampleData->Load_cal_RN_POI_BN();
				newSampleData->Load_Twitter_Social_Network_add_info(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true);

				GLOBAL_CHECKED_q.reset();
				GLOBAL_TESTED_POS_q.clear();
				GLOBAL_CHECHED_TESTED_q_COUNTER = 0;

				RunAllPossibleParameteres(newSampleData, Q_sets, filename, true, is_pruning_pow_eva, 1, load_tree);
				delete newSampleData;
			}
			catch (const std::exception& e)
			{
				Print_To_File(filename + "_error.txt", e.what());
			}

			break;
		case 2:
			try {
				filename = "Epinions";
				newSampleData = new SampleData(75888, 104770, no_p_r, no_p_s, _max_x, _max_y, 21048, total_no_of_keys);
				newSampleData->Load_cal_RN_POI_BN();
				newSampleData->Load_Epinions_Social_Network_add_info(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true);

				GLOBAL_CHECKED_q.reset();
				GLOBAL_TESTED_POS_q.clear();
				GLOBAL_CHECHED_TESTED_q_COUNTER = 0;
				RunAllPossibleParameteres(newSampleData, Q_sets, filename, true, is_pruning_pow_eva, 1, load_tree);
				delete newSampleData;
			}
			catch (const std::exception& e)
			{
				Print_To_File(filename + "_error.txt", e.what());
			}

			break;
		case 3:
			try {
				filename = "dblp_un_graph";
				newSampleData = new SampleData(317080, 104770, no_p_r, no_p_s, _max_x, _max_y, 21048, total_no_of_keys);
				newSampleData->Load_cal_RN_POI_BN();
				newSampleData->Load_dblp_ungraph_add_info(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true);
				GLOBAL_CHECKED_q.reset();
				GLOBAL_TESTED_POS_q.clear();
				GLOBAL_CHECHED_TESTED_q_COUNTER = 0;
				RunAllPossibleParameteres(newSampleData, Q_sets, filename, true, is_pruning_pow_eva, 1, load_tree);
				delete newSampleData;
			}
			catch (const std::exception& e)
			{
				Print_To_File(filename + "_error.txt", e.what());
			}
			break;
		default:
			std::cout << "Invalid Option" << std::endl;
			break;
		}
		break;

	case 2:
		try {
			int distribution_type = 0;
			std::cout << "distribution type ? Type: (1)-Uniform, (2)-Gaussian, (3)-Zipf " << std::endl;
			std::cin >> distribution_type;
			newSampleData = new SampleData(number_users[0], number_POI[default_number_POI_index], no_p_r, no_p_s, _max_x, _max_y, number_intersection_points[default_number_Intersection_points_index], total_no_of_keys);
			newSampleData->Generate_POI(max_no_of_keys_per_p, distribution_type);
			std::string distribution_type_name = distribution_type == 1 ? "uniform_" : distribution_type == 2 ? "gaussian_" : "zipf_";
			newSampleData->Print_poi(true, distribution_type_name + "ip_" + std::to_string(number_intersection_points[default_number_Intersection_points_index]) + "_poi_" + std::to_string(number_POI[default_number_POI_index]) + ".txt");
			delete newSampleData;
		}
		catch (const std::exception& e)
		{
			Print_To_File(filename + "_error.txt", e.what());
		}

		break;
	case 3:
		try {
			int _min_no_follower = 8;
			_max_no_follower = 40;
			int distribution_type = 0;
			std::cout << "distribution type ? Type: (1)-Uniform, (2)-Gaussian, (3)-Zipf " << std::endl;
			std::cin >> distribution_type;
			newSampleData = new SampleData(number_users[0], number_POI[default_number_POI_index], no_p_r, no_p_s, _max_x, _max_y, number_intersection_points[default_number_Intersection_points_index], total_no_of_keys);
			filename = "BiPartiteGraph";
			newSampleData->generate_Socail_Netword_with_k_d_truss(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true, filename, number_users, _min_no_follower, _max_no_follower, distribution_type);
			delete newSampleData;
		}
		catch (const std::exception& e)
		{
			Print_To_File(filename + "_error.txt", e.what());
		}

		break;
	case 4:
		try {

			filename = "BiPartiteGraph";

			int start_i = 0;

			int distribution_type = 0;
			std::cout << "distribution type ? Type: (1)-Uniform, (2)-Gaussian, (3)-Zipf , (4)-All " << std::endl;
			std::cin >> distribution_type;
			std::cout << "Start from ? Type: (0)-10k, (1)-20k, (2)-30k, (3)-40k (4)-50k (5)-100k (6)-200k" << std::endl;
			std::cin >> start_i;
			int stop = 0;
			std::cout << "Stop after this ? Type: (0)-10k, (1)-20k, (2)-30k, (3)-40k (4)-50k (5)-100k (6)-200k " << std::endl;
			std::cin >> stop;

			int start_j = distribution_type == 4 ? 1 : distribution_type;
			int stop_j = distribution_type == 4 ? distribution_type : distribution_type + 1;
			for (int i = start_i; i < 7; i++) {//vary No_of_users
				for (int j = start_j; j < stop_j; j++) {
					newSampleData = new SampleData(number_users[i], number_POI[default_number_POI_index], no_p_r, no_p_s, _max_x, _max_y, number_intersection_points[default_number_Intersection_points_index], total_no_of_keys);
					newSampleData->Load_BN(filename, false, j);

					if (i == start_i) {
						GLOBAL_CHECKED_q.reset();
						GLOBAL_TESTED_POS_q.clear();
						GLOBAL_CHECHED_TESTED_q_COUNTER = 0;
					}
					RunAllPossibleParameteres(newSampleData, Q_sets, filename, i != default_number_users_index, is_pruning_pow_eva, j, load_tree);
					delete newSampleData;
				}
				if (stop == i) break;
			}
		}
		catch (const std::exception& e)
		{
			Print_To_File(filename + "_error.txt", e.what());
		}
		break;
	case 5:
		try {


			std::vector<int> Q = Q_sets[default_Q];
			int k_value = k[default_k];
			int d_value = d[default_d];
			float omega_value = omega[default_omega];
			float pi_value = pi[default_pi];
			float theta_value = theta[default_theta];
			float sigma_value = sigma[default_sigma];
			int number_of_users = 30000;
			Q_sets[default_Q] = { 27, 47, 34, 20, 44, 13, 2 };
			float bs_w = 1, ss_w = 2, rs_w = 1;
			float stability_margin = 0.05;
			std::vector<int> user_updates = { 10,15,25,50,100,200 };

			vector<vector<int>> picked_users_for_each;
			for (int n : user_updates)
				picked_users_for_each.push_back(Pick_random_users_for_dynamic(n, number_of_users));

			int min_day = 50;
			int max_day = 60;


			total_no_of_keys = 50;
			filename = "BiPartiteGraph";

			DynamicApproch* newDynamicApproch = new DynamicApproch();
			std::unordered_set<int> comm_indecies_to_test= newDynamicApproch->Select_random_comm_index_to_test(filename,  5);


			std::ostringstream res;

			res << "user_updates,users_affcted,Repair_Dataset_Time,Offline_Recalculation_Time,Community_Load_Time,comm_count,Total_Community_Repair_Time,Avg_Community_Repair_Time,Total_Community_Repair_Time_batch,Avg_Community_Repair_Time_batch,Tree_Node_Repair_Time,Tree_Node_Repair_margin\n";
			
			for (int i = 0; i < user_updates.size(); i++) {
				newSampleData = new SampleData(number_of_users, number_POI[default_number_POI_index], no_p_r, no_p_s, _max_x, _max_y, number_intersection_points[default_number_Intersection_points_index], total_no_of_keys);
				newSampleData->Load_BN(filename, false, 1);				
				res << newDynamicApproch->Run_Dynamic_Temporal_Insert_visits(newSampleData, user_updates[i], picked_users_for_each[i], 5, min_day, max_day, filename, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value,bs_w,ss_w,rs_w, stability_margin, comm_indecies_to_test) << "\n";
								
				delete newSampleData;
			}
			Print_To_File_main(filename + "_dynamic_temporal_insert_percent_results.csv", res.str());


			res.clear();
			res << "user_updates,users_affcted,Repair_Dataset_Time,Offline_Recalculation_Time,Community_Load_Time,comm_count,Total_Community_Repair_Time,Avg_Community_Repair_Time,Total_Community_Repair_Time_batch,Avg_Community_Repair_Time_batch,Tree_Node_Repair_Time,Tree_Node_Repair_margin\n";
		
			for (int i = 0; i < user_updates.size(); i++) {
				newSampleData = new SampleData(number_of_users, number_POI[default_number_POI_index], no_p_r, no_p_s, _max_x, _max_y, number_intersection_points[default_number_Intersection_points_index], total_no_of_keys);
				newSampleData->Load_BN(filename, false, 1);
				
				res << newDynamicApproch->Run_Dynamic_Temporal_Delete_visits_by_user_number(newSampleData, user_updates[i], picked_users_for_each[i], 5, max_day, filename, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value, bs_w, ss_w, rs_w, stability_margin, comm_indecies_to_test) << "\n";

				delete newSampleData;
			}
			Print_To_File_main(filename + "_dynamic_temporal_deletion_percent_results.csv", res.str());
			delete newDynamicApproch;

			
			Q_sets[default_Q] = { 53,42,52,66,0,6,4 };
			total_no_of_keys = 72;
			filename = "Twitter";
			DynamicApproch* newDynamicApprochReal = new DynamicApproch();
			std::unordered_set<int> comm_indecies_to_test_real = newDynamicApprochReal->Select_random_comm_index_to_test(filename, 5);

			std::ostringstream resTw;

			resTw << "user_updates,users_affcted,Repair_Dataset_Time,Offline_Recalculation_Time,Community_Load_Time,comm_count,Total_Community_Repair_Time,Avg_Community_Repair_Time,Total_Community_Repair_Time_batch,Avg_Community_Repair_Time_batch,Tree_Node_Repair_Time,Tree_Node_Repair_margin\n";
			
			for (int i = 0; i < user_updates.size(); i++) {
				newSampleData = new SampleData(81306, 104770, no_p_r, no_p_s, _max_x, _max_y, 21048, total_no_of_keys);
				newSampleData->Load_cal_RN_POI_BN();
				newSampleData->Load_Twitter_Social_Network_add_info(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true);
				
				resTw << newDynamicApprochReal->Run_Dynamic_Temporal_Insert_visits(newSampleData, user_updates[i], picked_users_for_each[i], 5, min_day, max_day, filename, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value, bs_w, ss_w, rs_w, stability_margin, comm_indecies_to_test_real) << "\n";
				delete newSampleData;
			}
			Print_To_File_main(filename + "_dynamic_temporal_insert_results.csv", resTw.str());


			resTw.clear();
			res << "user_updates,users_affcted,Repair_Dataset_Time,Offline_Recalculation_Time,Community_Load_Time,comm_count,Total_Community_Repair_Time,Avg_Community_Repair_Time,Total_Community_Repair_Time_batch,Avg_Community_Repair_Time_batchime,Tree_Node_Repair_Time,Tree_Node_Repair_margin\n";
			
			for (int i = 0; i < user_updates.size(); i++) {
				newSampleData = new SampleData(81306, 104770, no_p_r, no_p_s, _max_x, _max_y, 21048, total_no_of_keys);
				newSampleData->Load_cal_RN_POI_BN();
				newSampleData->Load_Twitter_Social_Network_add_info(max_no_checkin_loc, max_euclidean_dist_btwn_user_POI, max_frequency, true);			
				resTw << newDynamicApprochReal->Run_Dynamic_Temporal_Delete_visits_by_user_number(newSampleData, user_updates[i], picked_users_for_each[i], 5, max_day, filename, Q, k_value, d_value, omega_value, pi_value, theta_value, sigma_value, bs_w, ss_w, rs_w, stability_margin, comm_indecies_to_test_real) << "\n";
				delete newSampleData;
			}
			Print_To_File_main(filename + "_dynamic_temporal_deletion_results.csv", resTw.str());
			delete newDynamicApprochReal;
		}
		catch (const std::exception& e)
		{
			Print_To_File(filename + "_error.txt", e.what());
			std::cout << "Error: " << e.what() << std::endl;
			int c;
			std::cin >> c;
		}
		break;
	default:
		std::cout << "Invalid Option" << std::endl;
		break;
	}
	return 0;

}