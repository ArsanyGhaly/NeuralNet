#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstddef>
#include <cmath>
#include <time.h>
#include <omp.h>
#include "NeuralNetArchitecture.h"
#include <numeric>
#include <chrono>

using namespace std::chrono;
using namespace std;
vector<string> seprate(const string &s, char sperator);


int main(int argc, char *argv[]){
	//enter the name of the file which cotain the data
	if (argc < 2){

		cout << "Error: Please enter the Name of the file" << endl;
		exit(EXIT_FAILURE);
	}

	//get the name and try to open 
	fstream inFile(argv[1]); 
	if (!inFile){
		cout << "Error:there is no such File " << argv[1] << endl;
		exit(EXIT_FAILURE);

	}

	//if there is no error get data file name
	string temp;
	getline(inFile, temp, '\n');
	cout <<"the name of the file: " << temp <<endl;
	
    //get the number of the hidden layers
	unsigned int hidden;
	inFile >> hidden;
	int Layernodes[hidden + 2];
	cout <<"the hidden layers : "<< hidden <<endl;

	//get number of nodes in each hidden layer 
	for (int i = 1; i <= hidden; ++i){
		inFile >> Layernodes[i];
		cout <<"layer : "<<i <<"nodes" <<Layernodes[i]<< endl;
	}

	//get the learning rate and error tolerance
	double learningRate;
	double errorTolerance;
	inFile >> learningRate;
	inFile >> errorTolerance;
	cout <<"learningRate : "<< learningRate <<endl;
	cout <<"errorTolerance : "<< errorTolerance <<endl;

	//done with is file
	inFile.close();

	//open the data file and start reading
	fstream dataFile(temp.c_str());

	if (!dataFile){
		cout << "Error: Can not open such a data file: " << temp << endl;
		exit(EXIT_FAILURE);

	}

	
	//now get the data from the datafile
	double scaler;
	dataFile >> scaler;
	cout <<"scaler : "<< scaler <<endl;

	//get the number of input and the number of out put
	string io;
	dataFile >> io;
	
	vector<string> inputs_string;
	int count =0;
	while(dataFile >> temp){
		count++;
		if(temp[0] == 0) continue;
		inputs_string.push_back(temp);
	
	}
	int num_nodes = count;
	cout <<"number of Nodes: "<< num_nodes <<endl;
	dataFile.close(); 

	vector<string> strtemp = seprate(io, ',');
    unsigned int i_num = stoul(strtemp[0]);
    unsigned int o_num = stoul(strtemp[1]);
    Layernodes[0] = i_num;
    Layernodes[hidden+1] =o_num;
   	
   	cout <<"number of input : "<< i_num <<endl;
   	cout <<"number of output : "<< o_num <<endl;
    
  
	vector<double*> data;

	for (int i = 0; i < num_nodes; ++i){
		vector<string> temp = seprate(inputs_string[i], ',');

		double *temp2 = new double[i_num + o_num];

		for (int j = 0; j < i_num + o_num; ++j){

			temp2[j] = atof(temp[j].c_str());
			

		}

		data.push_back(temp2);
		cout <<data.size();

	}

NeuralNetArchitecture net(errorTolerance, learningRate, scaler,i_num,o_num , hidden, Layernodes,num_nodes);    
net.optamization(data,10,num_nodes);




/*
	//to keep track how many loop
  
    double *answer[num_nodes];

    //inizlize the answer array
	for (int i = 0; i < num_nodes; ++i)
	{
		answer[i] = new double[1];
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	

	//push the answers
	  int run = 0;
	do{
			for (int i = 0; i < num_nodes; ++i){

				delete [] answer[i];
				answer[i] = net.pro_back(data[i]);

			}
			run++;
			cout <<run;

	}while((duration_cast<duration<double>>(high_resolution_clock::now() - t1)).count() <= errorTolerance);

cout<<"=================================================="<<endl;
cout <<"running for:"<<run <<endl;
cout << "For each Point:"<<endl;
double error = 0.0;
	vector<vector<double> > actual_out; 
	vector<vector<double> > our_out;

for (int i = 0; i < num_nodes; ++i){
		vector<double> our; 
		vector<double> actual;
	    cout << "for Node" << i + 1 << ": " << endl;
		for (int j = 0; j < i_num; ++j){

				cout << " " << data[i][j]<<endl;

			}

			cout << endl << "Our outputs:  ";

			for (int j = 0; j < o_num; ++j){

				cout << answer[i][j]<<endl;
				our.push_back(answer[i][j]);

			}

			cout << endl << "Real outputs:";

			for (int j = 0; j < o_num; ++j){

				cout << data[i][i_num + j]<<endl;
				actual.push_back(data[i][i_num + j]);

			}

			actual_out.push_back(actual);
			our_out.push_back(our);

			

			cout << "================================"<< endl;
		}
		
	for (int i = 0; i < num_nodes; ++i){
		for(unsigned int j = 0; j < o_num; j++)
			{
			double e = abs(actual_out[i][j] - (our_out[i][j] * scaler)) / actual_out[i][j];
			error += e;
			}
		}
		//free the answer array for the next round
		for (int i = 0; i < o_num; ++i){

			delete [] answer[i];

		}
 	cout << "Error: " << (double)(error * 100) << "%" <<endl;
	cout << "Avg error: " << (double)(error * 100 / num_nodes) << "%" <<endl;
	*/
	return 0;
}


//total

		
	

















vector<string> seprate(const string &s, char sperator){
	stringstream temp(s);
	string item;
	vector<string> items;
	while(getline(temp, item, sperator)){
		items.push_back(item);
	}
	return items;
}