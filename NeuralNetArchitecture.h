#ifndef NEURAL_NET_ARCHITECTURE_H
#define NEURAL_NET_ARCHITECTURE_H

#include <stdexcept>
#include <exception>
#include <cmath>
#include <iostream>
#include <vector>


#include "NeuralNet.h"
using namespace std;

class NeuralNetArchitecture {
	
private:
	//layers
	vector<NeuralNet*> N_nets;
	int *layersize;				 //eachlayer size
	unsigned int hiddenlayers;   //number of the hidden layers

	//helpers
	double Tolerance;
	double learningRate;         
	double scaler;               //scaler

	//sizes
	unsigned int input_size;    //1
	unsigned int output_size;   //1
	int solve_point(const double x);
	void train(const vector<double*> &data_inputs, NeuralNet *n);
public:
	NeuralNetArchitecture(double T, double R, double M,unsigned int in, unsigned int out, unsigned int H, int *Ls,unsigned int nodes_num);
	~NeuralNetArchitecture();
	void optamization(vector<double*> data,unsigned int bins_size,unsigned int nodes_num);
	vector<vector<vector<double> > >  k_fold_cross(vector<double*> data,unsigned int bins_size);
	double** trainNeuralNet(vector<vector<double> > data, NeuralNet *net, unsigned int num_nodes);
	double test_network(NeuralNet* net,vector<vector<double> > data);
	vector<vector< int> > get_nodes();

};

#endif
