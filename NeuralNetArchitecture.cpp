#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstddef>
#include <cmath>
#include <time.h>
#include <omp.h>
#include "NeuralNet.h"
#include <numeric>
#include <chrono>
#include <random>


#include "NeuralNetArchitecture.h"
using namespace std;
using namespace std::chrono;



//get Random numbers for set the wights up
double RandomInt(const int &x,const int &y){

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   default_random_engine generator(seed);
   std::uniform_int_distribution<unsigned int> distribution(x,y);
   
   return distribution(generator);

}
NeuralNetArchitecture::NeuralNetArchitecture(double T, double R, double M,unsigned int in, unsigned int out, unsigned int H, int *Ls,unsigned int nodes_num)
{
	//set instance veriables
	input_size = in;
	output_size =out;
	
	Tolerance = T;
	learningRate= R;
	scaler=M;

	//layersize = Ls;
	//hiddenlayers=H;
	//for(int i=0; i < nodes_num; i++){
		//unsigned int size = points[i].size();
		//int* p =(int*) points[i].data();
		//cout <<"done1" <<endl;
	   // NeuralNet* temp = new NeuralNet(Tolerance, learningRate, scaler,input_size, output_size,hiddenlayers, layersize);
	    //N_nets.push_back(temp);
	//}

}
	
NeuralNetArchitecture::~NeuralNetArchitecture(){
}

double** NeuralNetArchitecture::trainNeuralNet(vector<vector<double> > data, NeuralNet *net, unsigned int num_nodes){
 	double *answer[num_nodes];


    //inizlize the answer array
	for (int i = 0; i < num_nodes; ++i)
	{
		answer[i] = new double[1];
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
	double temp[num_nodes][input_size + output_size];
		for (int i = 0; i < num_nodes; ++i){
			for (int j = 0; j < input_size + output_size; ++j){

						temp[i][j] =data[i][j];
						

					}


				}
	//push the answers
	  int run = 0;
	do{
			for (int i = 0; i < num_nodes; ++i){

				delete [] answer[i];
				answer[i] = net->pro_back(temp[i]);

			}
			run++;
			//cout <<run;

	}while((duration_cast<duration<double>>(high_resolution_clock::now() - t1)).count() <= Tolerance);
return answer;
}


double NeuralNetArchitecture::test_network(NeuralNet* net,vector<vector<double> > data)
{
	int data_size = data.size();
	double **answer =trainNeuralNet(data,net,data_size);

    vector<vector<double> > actual_out; 
	vector<vector<double> > our_out;
	double error = 0.0;
/*
	for(unsigned int i = 0; i < data_size; i++){
	    vector<double> our; 
		vector<double> actual;
		
		//our output
		for(unsigned int j = 0; j < output_size; ++j){ 
			our.push_back(answer[i][j]);
		}

		for(unsigned int j = 0; j < output_size; j++){ 
			actual.push_back(data[i][j+ input_size]);
		}
	actual_out.push_back(actual);
	our_out.push_back(our);
	}

	for(unsigned int i = 0; i < data_size; i++){
		
		for(unsigned int j = 0; j < output_size; j++){
			double temp = abs(actual_out[i][j] - (our_out[i][j] * scaler));
			error += (temp/actual_out[i][j]);
		}
	}*/
	return error;
}

vector<vector<vector<double> > >  NeuralNetArchitecture::k_fold_cross(vector<double*> data,unsigned int bins_size)
{
//bins >>eachbin>>eachvector on the bin
vector<vector<vector<double> > > bins;

//inislize the bins
for(unsigned int i = 0; i < bins_size; i++){
		vector<vector<double> > inner;
		bins.push_back(inner);
	}
//start spliting the data
	unsigned int counter = 0; 
	vector<unsigned int> visited;

// loop through each unvisited data.
while(visited.size() < data.size()){
		unsigned int index = (unsigned ) RandomInt(0, (data.size() - 1));
		//go over the visited data set 
		bool f= false;
		//if the index is visited get new rondom index
		for(unsigned int j = 0; j < visited.size(); j++){
			if(visited[counter] == index){ 
				f = true;
				break;
			}
		}
		//if Not containue
		if(f){
			continue;
		}

		visited.push_back(index);

		// hold the input and the output
		vector<double> IO_hold;
		for(unsigned int i = 0; i < (input_size + output_size); i++){
			IO_hold.push_back(data[index][i]);
		}


		//put the input/out put factor in a bin.
		bins[counter].push_back(IO_hold);
		cout <<bins[counter].size();
		
		// Increment bin_counter, if it e = num_bins, ==> 0.
		counter++;
		counter %= bins_size;
		//cout<< "round"<<counter <<endl;
	}
	
	return bins;

}


vector<vector<int> > NeuralNetArchitecture::get_nodes(){
	vector<vector< int> > t;
	for(int i = 0; i < 5; i++){
		vector< int> temp;
		
		for(int j = 0; j <i; j++){
			temp.push_back(5);
		}
		
	
		int total = pow(5, i);
		
		
		for(int k = 0; k < total; k++){
			vector< int> tt(i+2);
			
			tt[0] = input_size;
			tt[i+1] = output_size;
			
			for(int j = 1; j <= i; j++){
				tt[j] = temp[j-1];
			}
			
			// If we are not on layer 0
			bool flag = false;
			
			if(i > 0){
				
				// Decrement a node in a layer
				int in = i-1;
				temp[in]--;
				
				while(temp[in] == 0){
					
					if(in == 0 && temp[in] == 0){
						flag = true;
						break;
					}
					temp[in] = 10;
					in--;
					temp[in]--;
				}
			}
			t.push_back(tt);
			if(flag) break;
		}
	}
	return t;
}

void NeuralNetArchitecture::optamization(vector<double*> data,unsigned int bins_size,unsigned int nodes_num)
{
	vector<vector<vector<double> > > bins = k_fold_cross(data,bins_size);
	
	vector<vector<int> > points = get_nodes();
	
    int** lss;
	cout << "points.size()"<<points.size()<<endl;
	for (int i =0;i<points.size();i++){
		 int *ls_1;
		//cout << "folds"<<endl;
		for (int j =0;j<points[i].size();j++){
		ls_1[j] = points[i][j];
		
	}

	lss[i] = ls_1;
		}


	for(int i =0; i< 30;i++){
		int size = points[i].size();
		NeuralNet* temp = new NeuralNet(Tolerance, learningRate, scaler,input_size, output_size,size,lss[i]);
	    N_nets.push_back(temp);
		cout <<"=================="<<endl;
	}
	

	//vector<NeuralNet*> networks;
	


	

	double init_error = 99999.9;
    int nn =1;
    int total=0;
    unsigned int *ls;
    int net_id =0;

for(unsigned int c=0;c<N_nets.size()-1;c++)
{
	
	double error_cur = 0.0;
	//loop through all the bins
	double t[nodes_num][input_size + output_size];
	for(unsigned int i = 0; i < bins_size; i++){
			int rand_bin =RandomInt(0, bins_size-1);
			vector<vector<double> > inputs;
			vector<vector<double> > outputs;
		//choose a rondom bi
			
			for(unsigned int j = 0 ; j < bins_size; j++){
			//if we choose it before containue
			if(j == rand_bin) continue;
			
			//int total = bins[j].size();
			
			for(unsigned int k = 0; k < bins[j].size(); k++){
				
				//hold the input and the output for this bin
				vector<double> in_data(input_size);
				vector<double> out_data(output_size);
				
				//push the inputs for that bin 
				for(unsigned int h  = 0; h < input_size; h++)
				{
					in_data[h] = bins[j][k][h];
					t[k][h]=bins[j][k][h];
					
				}
				
				//push the outpit for that bin		
				for(unsigned int o = 0; o < output_size; o++)
				{
					out_data[o] = bins[j][k][o+input_size];
					t[k][input_size + o]=bins[j][k][o+input_size];
					
				}
					//get the in put and output to construct 
				    inputs.push_back(in_data);
				    outputs.push_back(out_data);
					
				}
			}
			    
		
				cout <<"done"<<endl;
				high_resolution_clock::time_point t1 = high_resolution_clock::now();
	

			    N_nets[c]->pro_back(t[c]);
				test_network(N_nets[c], bins[rand_bin]);
			




		
				

				}

				//train this net with the data of this inout and output
				
			
			}

			//calculate the error for each net
			//error_cur = (error_cur / bins_size) * 100;
cout <<"done100"<<endl;		

}








	







