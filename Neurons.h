#ifndef NEURONS_H_
#define NEURONS_H_
#include<vector>
#include <iostream>

using namespace std;
class Neurons {
		private:
			int id;
			double value;                                   
			double wightSum;	                              
			double delta;		
			double input_bias;							
			
		public:
			//constructors
				Neurons(unsigned int id);
				~Neurons();
				
			//setters
				bool set_id(int id);
				bool set_bias(double bias);
				bool set_value(double value);
				bool set_delta(double d);
				bool set_sum(double s);
				
			
			//getters
				int get_id();
				double get_bias();
				double get_delta();
				double get_value();
				double get_sum();

			//print the node 
				void prient_node(); 
				//get Random numbers for set the wights up



};





#endif 