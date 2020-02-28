//Author: Thomas Merkh, RPI Class of 2016, merkht@rpi.edu//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream>
#include "rng.cpp"

using namespace std;

typedef long double ld;

int main()
{

	//seed random number generators (dont touch)
	srand(time(NULL));
	RandomInitialise(rand()%30000, rand()%30000);
	
	//initialize constants
	ld pi = 3.14159265359;
	ld K = 0.46;                 //
	ld lambda = 39;              //nm
	ld w = 40;                   //width
	ld h = 25;                   //height
	ld A_ = h*w;                 //mean area
	ld Constant = (1  +  K*((lambda/w) + (lambda/h)));

	

	//dont need to touch
	const int N = 10000000;      //number of points to use in integration
	const int ints = 5;          //number of integrations to average to get final resistivity (higher integer reduces numerical error)
	ld single_summ = 0;			 //integral value without averaging
	ld final_summ = 0;			 //final integral value to be averaged
	ld upper_bound = 100000;     //upper limit --> DONT INCREASE due to numerical error
	ld ep = 0.01;                // 0 < ep << 1 (lower integration limit)
	ld A, rho_rho0, constant1;              //Initializing the variables
	ofstream out("Line_edge.csv");


	for (ld s = 1; s <= 50; ++s)
	{
		constant1 = (A_) / (s*sqrt(2*pi));
		rho_rho0 = 0;
		final_summ = 0;
		for (int j = 1; j <= ints; ++j)
		{
			single_summ = 0;
			for (int i = 0; i < N; ++i)
			{
				A = RandomDouble(ep,upper_bound);
				single_summ += (1/A)*exp( -1*((A - A_)*(A - A_)) / (2*s*s) );
			}
			final_summ += (upper_bound-ep)*single_summ/N;
			cout << j << " integrations preformed for s = " << s << endl; 
		}
		//final integral value
		rho_rho0 = final_summ/ints;
		rho_rho0 = rho_rho0*Constant*constant1;
		cout << "Resistivity for s = " << s << ", is: " << rho_rho0 << endl;
		out << s << "," << Constant << "," << constant1 << "," << final_summ/ints << ","  << rho_rho0 << endl;
	}
	
    out.close();
	return 0;
}