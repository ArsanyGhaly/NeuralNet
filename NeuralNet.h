#ifndef NEURAL_NET_H
#define NEURAL_NET_H

#include <stdexcept>
#include <exception>
#include <cmath>
#include <iostream>
#include <vector>

#include "Neurons.h"

using namespace std;

class NeuralNet {
	
private:
	//layers
	Neurons ***layers;           //layers contains nodes
	vector<Neurons*> neurons;    //vector of the whole nodes
	int *layersize;				 //eachlayer size
	unsigned int hiddenlayers;   //number of the hidden layers

	//helpers
	double Tolerance;
	double learningRate;         
	double scaler;               //scaler

	//sizes
	unsigned int input_size;    //1
	unsigned int output_size;   //1
	unsigned int totle;       

	//wights for nodes
	double *segments_wights;  

	//tools
	
	double sigmode(double in);
	double compute_delta(double in );




public:

	
	NeuralNet(double T, double R, double M,unsigned int in, unsigned int out, unsigned int H, int *Ls);
	
	
	~NeuralNet();
	double *Pro_forward(double *z);
	double *pro_back(double *z);
	
};

#endif
