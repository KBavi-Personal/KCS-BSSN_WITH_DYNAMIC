//--------------------------------------------
// this function contains necessary random number
// generators
// 
// collected by Yufei Tao
//--------------------------------------------

#ifndef RAND_H
#define RAND_H
namespace random1 {
	//--------------------------------------------
	double gaussian_org(double mean, double sigma);
	double gaussian_pos(double min, double max, double mean, double sigma);
	double uniform_Yufei(double _min, double _max);
	double zipf_org(double x1, double x2, double p);
	double zipf_pos(double x1, double x2, double p);

	int uniform_int_org(int _min, int _max);
	int uniform_int_pos(int _min, int _max);
	//--------------------------------------------
}


#endif