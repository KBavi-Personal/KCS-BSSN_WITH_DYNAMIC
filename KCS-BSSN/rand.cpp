//--------------------------------------------
// this file contains necessary random number
// generators
// 
// collected by Yufei Tao
//--------------------------------------------
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "rand.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <random>
#include <chrono>

using namespace std;
/************************************************************
***  Given a mean and a standard deviation, gaussian       **
**   generates a normally distributed random number        **
**   Algorithm:  Polar Method, p.  104, Knuth, vol. 2      **
************************************************************/

double random1::gaussian_org(double mean, double sigma)
{
	double v1, v2;
	double s;
	double x;

	do
	{
		v1 = 2 * random1::uniform_Yufei(0, 1) - 1;
		v2 = 2 * random1::uniform_Yufei(0, 1) - 1;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.);

	x = v1 * sqrt(-2. * log(s) / s);

	/*  x is normally distributed with mean 0 and sigma 1.  */
	x = x * sigma + mean;

	return (x);
}
double random1::gaussian_pos(double min,double max,double mean, double sigma)
{
	double x=gaussian_org(mean, sigma);
	while (x < min || x > max)
		x = gaussian_org(mean, sigma);

	return (x);
}
/************************************************************
** Generates a random number between _min and _max         **
** uniformly                                               **
   By Yufei Tao
************************************************************/

double random1::uniform_Yufei(double _min, double _max)
{
	//	cout<<_min<<"  "<<_max<<endl;
	int int_r = rand();
	long base = RAND_MAX - 1;
	double f_r = ((double)int_r) / base;
	double rlt = (_max - _min) * f_r + _min;

	if (rlt > _max)
		rlt = _max;
	if (rlt < _min)
		rlt = _min;
	/*
		while (rlt<_min ||rlt>_max)
		{
	//		cout<<"rand error:"<<endl;

			int_r = rand();
			f_r  = ((double) int_r) / base;
			rlt =(_max - _min) * f_r + _min;

		}
	*/
	return rlt;
}

/*************************************************************/
/*  zipf generates a random number that follows Zipf         **
**  distribution and lies between x1 and x2.                 **
**  original code by Christos Faloutsos, 1995

**  The original node outputs discrete data only. The current**
**  function remedies the problem.			                 **
**  Modified by Yufei Tao (08/Dec/02)                         **
**************************************************************/
double random1::zipf_org(double x1, double x2, double p)
{

	double x;
	double i;
	double r, HsubV, sum;
	int V = 100;

	//double uniform();

	/* calculate the V-th harmonic number HsubV. WARNING: V>1 */
	HsubV = 0.0;
	for (i = 1; i <= V; i++)
		HsubV += 1.0 / pow((double)i, p);

	r = uniform_Yufei(0., 1.) * HsubV;
	sum = 1.0; i = uniform_Yufei(0, 1);
	while (sum < r) {
		//i++;  //commented by Yufei Tao
		i += uniform_Yufei(1, 2);
		sum += 1.0 / pow((double)i, p);
	}

	/* i follows Zipf distribution and lies between 1 and V */

	/* x lies between 0. and 1. and then between x1 and x2 */
	x = ((double)i - 1.) / ((double)V - 1.);
	x = (x2 - x1) * x + x1;

	return(x);
}

double random1::zipf_pos(double x1, double x2, double p)
{
	double x = zipf_org( x1,  x2,  p);
	while (x < x1 || x > x2)
		x = zipf_org(x1, x2, p);

	return (x);
}

int random1::uniform_int_org(int _min, int _max)
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 gen(seed);
	std::uniform_int_distribution<int> int_dist(_min, _max);
	return int_dist(gen);
}
int random1::uniform_int_pos(int _min, int _max)
{
	double x = uniform_int_org(_min, _max);
	while (x < _min || x > _max)
		x = uniform_int_org(_min, _max);

	return (x);
}

