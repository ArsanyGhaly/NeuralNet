#include "Neurons.h"
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include <chrono>
using namespace std;
double RandomD(const double &x,const double &y){

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   default_random_engine generator(seed);
   uniform_real_distribution<double> distribution(x, y);
   
   return distribution(generator);

}



Neurons::Neurons(unsigned int _id)
{
		id = _id;
		input_bias = RandomD(-0.1,0.1);
		delta = 0.0;
		wightSum = 0.0;
		value =0.0;
		
}

Neurons::~Neurons()
{}


bool Neurons::set_id(int _id)
{
	id=_id;
}


bool Neurons::set_bias(double bias)
{
	input_bias = bias;
	return true;
}

bool Neurons::set_value(double v)
{
	value = v;
	return true;
}

bool Neurons::set_delta(double d)
{
	delta = d;
	return true;
}
bool Neurons::set_sum(double s){
	wightSum = s;
}

double Neurons::get_delta()
	{
		return delta;
	}

double Neurons::get_bias()
	{
		return input_bias;
	}
 
double Neurons::get_value()
	{
		return value;
	}

int Neurons::get_id()
	{
		return id;
	}


double Neurons::get_sum(){
	return wightSum;
}



void Neurons::prient_node(){
	cout<<get_id() <<get_bias() << get_delta() <<get_value() <<get_sum() <<endl;
} 



