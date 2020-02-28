//Author: Thomas Merkh
//Rensselaer Mathematics Department, Class of 2016

#ifndef Particle_h_
#define Particle_h_

#include "rng.cpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <utility>

 //<cstdlib>

using namespace std;

class Particle
{
public:
	Particle(){};
	Particle(int i, double ymax_){
		number = i;
		x = 0;
		y = 0;
		ymax = ymax_;
		steps = 0;
	};

	//public functions
	void move(double lambda);
	double deposit(double lambda, int a, int b);

private:
	int number;
	double x;
	double y;
	double ymax;
	int steps;	                   
};

double Particle::deposit(double lambda, int a, int b){
	RandomInitialise(a%30000,b%30000);
	while(y < ymax){
		move(lambda);
	}
	//cout << "Particle number " << number << " took " <<  steps << " steps to deposit " << endl;
	return x;
}

void Particle::move(double lambda){
	double theta = RandomDouble(0,3.14);
	//double theta = RandomDouble(-3.14*5/12.0,3.14*17/12.0);  unbiased version
	double lamb = lambda + RandomGaussian(0,lambda/3.0);  //slightly more realistic mfp  //set arg to (0,0) for no sigma around lambda
	x = x + lamb*cos(theta);
	y = y + lamb*sin(theta);

	//if (steps%1000000 == 0){cout << y-ymax << endl;}  // moving backward
	steps++;
}

#endif
