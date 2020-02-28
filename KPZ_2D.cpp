#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "rng.cpp"

using namespace std;

int main()
{
	srand (time(NULL));
	RandomInitialise(rand()%30000,rand()%30000);
	vector<vector<double> > h;
	const int Nmax = 512;
	double v = 1;
	double lambda = -0.5;
	double mean = 0.0;  
	double stddev = 0.75;
	const int tmax = 10000;
	double L = 128;
	double dd = L/Nmax;
	double dt = 10.0/tmax;

	for (int i = 0; i < Nmax; ++i)
	{
		h.push_back(vector<double>());
		for (int j = 0; j < Nmax; ++j)
		{
			h[i].push_back(0);
		}
	}

	//set boundary conditions
	//for (int i = 0; i < Nmax; ++i)
	//{h[i][0] = 0;}

	int x,xx,y,yy;
	double avg, rms;
	ofstream roughness("RMS.csv");

	for (int t = 0; t < tmax-1; ++t)
	{
		if(t % (tmax/10) == 0){cout << "Time step: " << t << endl;}
		avg = 0;
		rms = 0;
		//update the surface heights
		for (int i = 0; i < Nmax; ++i)
		{
			for (int j = 0; j < Nmax; ++j)
			{
				//periodic conditions
				x = i+1;
				if(x == Nmax){x = x - Nmax;}
				xx = i-1;
				if(xx == -1){xx = xx + Nmax;}
				y = j+1;
				if(y == Nmax){y = y - Nmax;}
				yy = j-1;
				if(yy == -1){yy = yy + Nmax;}

				h[i][j] = h[i][j] + dt*( v*(h[x][j] + h[xx][j] + h[i][y] + h[i][yy] - 4*h[i][j])/(2*dd*dd) + (lambda/2)*(pow((h[x][j] - h[xx][j] + h[i][y] - h[i][yy])/(4*dd),2)) + RandomGaussian(mean,stddev));
			}
		}

		for (int i = 0; i < Nmax; ++i)
		{
			for (int j = 0; j < Nmax; ++j)
			{
				avg += h[i][j];
			}
		}
		avg = avg/(Nmax*Nmax);  //mean surface height

		for (int i = 0; i < Nmax; ++i)
		{
			for (int j = 0; j < Nmax; ++j)
			{
				rms += pow(avg - h[i][j],2);
			}
		}
		rms = sqrt((rms/(Nmax*Nmax)));

		roughness << t << "," << rms << endl;
	}
	roughness.close();

	ofstream output("KPZ.csv");
	for (int i = 0; i < Nmax; ++i)
	{
		for (int j = 0; j < Nmax; ++j)
		{
			output << h[i][j] << ",";
		}
		output << endl;
	}
	output.close();


	return 0;
}