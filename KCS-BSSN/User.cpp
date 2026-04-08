#include "User.h"
User::User(int _id, int _P_s_size, int total_n_keys)
{
	ub_sub = 0;
	ub_w_in = 0;
	ub_w_out = 0;
	id = _id;
	P_s_size = _P_s_size;
	dist_P_s.resize(P_s_size, -1);

	// for pruning
	isPruned = false;

	Keys_Fsum.resize(total_n_keys, 0.0);
	Keys_Fmax.resize(total_n_keys, 0.0);
	Keys_visited_count.resize(total_n_keys, 0);

	ub_f_sum = 0;
	ub_f_avg = 0.0f;
}
User::~User()
{
}


