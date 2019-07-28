#include <stdexcept>
#include <exception>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

#include "Neurons.h"
#include "NeuralNet.h"
using namespace std;


int Matrix_index(int r,int col,int size){
	int temp = (r *size) +col;
	return temp;
}

double RandomDouble(const double &x,const double &y){

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   default_random_engine generator(seed);
   uniform_real_distribution<double> distribution(x, y);
   
   return distribution(generator);

}



NeuralNet::NeuralNet(double T, double R, double M,unsigned int in, unsigned int out, unsigned int H, int *Ls)
{
	//set instance veriables
	input_size = in;
	output_size =out;
	totle =0;
	
	Tolerance = T;
	learningRate= R;
	scaler=M;

	layersize = Ls;
	hiddenlayers=H;
	
	cout <<"-----------------------------------------"<<endl;
 	cout <<"number of input : "<< input_size <<endl;
   	cout <<"number of output : "<< output_size <<endl;
   	cout <<"learningRate : "<< learningRate <<endl;
	cout <<"errorTolerance : "<< Tolerance <<endl;
	cout <<"scaler : "<< scaler <<endl;
	
	cout <<"the hidden layers : "<< hiddenlayers <<endl;
	//layer size input,hiddenlayers,output

	//put the Nodes in each layer
	//the layers will be the hidden layers + the input and the out
	layers = new Neurons**[hiddenlayers + 2];

cout <<"the hidden layers : "<< hiddenlayers <<endl;
	//inizlize the layers 
	for (int i = 0; i < hiddenlayers+2; ++i)
	{
		//enters the new layer
		cout <<"the hidden layers : "<< layersize[i] <<endl;
		layers[i] = new Neurons*[layersize[i]];
		
		//fill the layer with nodes
		for (int j = 0; j < layersize[i]; ++j){
			layers[i][j] = new Neurons(totle++);
			//layers[i][j]->set_bias(RandomDouble(-0.1, 0.1));
			//put the nodes inside the nodes vector
			neurons.push_back(layers[i][j]);

		}
		
	}
	


	//creat a matrix to store the wights of each node to each node wight
	//store 2d in 1d (use the (row * size) +col)
	segments_wights = new double[totle * totle];
	
	//get rondom wight for each node between -1 and 1
	for (int i = 0; i < totle; ++i)
		{		
				for(unsigned int j=0;j < totle; j++)
				{
					segments_wights[Matrix_index(i,j,totle)] = segments_wights[Matrix_index(j,i,totle)] = RandomDouble(-0.1, 0.1);

				}

				//change the wights to the new wights
				//neurons[i]->set_weights(new_weights, nwights); 
				//neurons[i]->set_bias(RandomDouble(-1, 1));
		}	

}
	
NeuralNet::~NeuralNet(){
for (int i = 0; i < hiddenlayers + 2; ++i){

	//free the node in each layer 
	for (int j = 0; j < layersize[i]; ++j)
	{
		delete layers[i][j];
	}
	//delete the layer itself
		delete [] layers[i];

	}
	//delete the layers
	delete []segments_wights;
	delete layers;

}

double *NeuralNet::Pro_forward(double *z)
{
//set the first layer to the inputs	
for (int i = 0; i < input_size; ++i)
	{
		layers[0][i]->set_value(z[i]);
		
		layers[0][i]->set_sum(z[i]);

	}
	

	//for every hidden layer but not the first
	for (int i = 1; i <= hiddenlayers + 1; ++i)
	{
		//for every sigle point in each layer
		for (int j = 0; j < layersize[i]; ++j)
		{
			Neurons *cur = layers[i][j];
			cur->set_sum(cur->get_bias());
			for (int h = 0; h < layersize[i - 1]; ++h)
			{
				Neurons *pre = layers[i - 1][h];
				double temp = pre->get_value() * segments_wights[Matrix_index(pre->get_id(),cur->get_id(),totle)];
				cur->set_sum(cur->get_sum() + temp);   
			}
		
		cur->set_value(sigmode(cur->get_sum()));
		}
	}
		
	double *out = new double[output_size];
	for (int i = 0; i < layersize[hiddenlayers + 1]; ++i)
	{

		out[i] = (layers[hiddenlayers + 1][i]->get_value() * scaler);

	}

	return out;
}

double *NeuralNet::pro_back(double *z)
{

	double *forword = Pro_forward(z);
	

	for (int i = 0; i < layersize[hiddenlayers + 1]; ++i)
	{
		Neurons *cur = layers[hiddenlayers + 1][i];
		double temp = (z[input_size + i] / scaler) - cur->get_value();
		cur->set_delta(compute_delta(cur->get_sum()) * temp);


	}


	//start at the end and go to each node in each layer 
	for(int l = hiddenlayers; l >= 0; --l)
	{
		for (int j = 0; j < layersize[l]; ++j)
		{
		Neurons *cur = layers[l][j];
		
		double summ = 0.0;
		
		for (int h = 0; h < layersize[l + 1]; ++h)
			{
			Neurons *n = layers[l + 1][h];
			double t = segments_wights[Matrix_index(cur->get_id(),n->get_id(),totle)];
			summ += (t * n->get_delta());

			}
			cur->set_delta(compute_delta(cur->get_sum())* summ);
		}
	}

	//update the wights the same way as we set them up 2d array
	for (int i = 0; i < totle; ++i){
		for (int j = 0; j < totle; ++j){
			segments_wights[Matrix_index(i,j,totle)]=segments_wights[Matrix_index(j,i,totle)]=segments_wights[Matrix_index(i,j,totle)] +(learningRate * neurons[j]->get_delta() * neurons[i]->get_value());;
			
	}
	neurons[i]->set_bias(neurons[i]->get_bias() + (learningRate * neurons[i]->get_delta()));
	
}
cout <<"back" <<endl;
return forword;

}


double NeuralNet::sigmode(double x){
	double tem = (1 / (1+exp(-x)));
	return tem;

}

double NeuralNet::compute_delta(double y ){
	double sig = sigmode(y);
	//dravative of the sigmoid
	double temp = (sig * (1 - sig));
	return temp;

}






