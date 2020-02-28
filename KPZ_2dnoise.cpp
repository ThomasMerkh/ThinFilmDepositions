//2-D improved noise

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <time.h>
#include <algorithm>
#include "rng.cpp"

using namespace std;
/* comment space...

*/


int main()
{
	srand (time(NULL));
	RandomInitialise(rand()%30000,rand()%30000);
	vector<vector<double> > h;
	vector<vector<double> > h_temp;
	const int L = 512;
	const double tmax = 5000;
	double dx = 1;
	double dt = 0.01; 
	double v = 1;
	double lambda = 10;
	double D = 1;
	double sigma = sqrt(2*D/dx);
	double g = lambda*lambda*D/(v*v*v);
	int x, xx, y, yy;
	double avg = 0;
	double rms = 0;
	double epsilon = 5*pow(10,-3);
	double one,two,three,four,five;

	//initialize size and initial condition
	for (int i = 0; i < L; ++i)
	{
		h.push_back(vector<double> (L,RandomDouble(0,1)));
		h_temp.push_back(vector<double>(L,RandomDouble(0,1)));  //L ints with value 0
	}

	ofstream rough("roughness.csv");
	for (int t = 0; t < tmax; ++t)
	{
		if (t%100 == 0){cout << "time: " << t << endl;}
		avg = 0;
		rms = 0;

		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < L; ++j)
			{
				//periodic conditions
				x = i+1;
				if(x == L){x = x - L;}
				xx = i-1;
				if(xx == -1){xx = xx + L;}
				y = j+1;
				if(y == L){y = y - L;}
				yy = j-1;
				if(yy == -1){yy = yy + L;}
				//interesting if it gets averaged over 3 occurences of noise?
				one = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				two = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				three = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				four = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				five = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				h_temp[i][j] = (one+two+three+four+five)/5.0;

				//h_temp[i][j] = h[i][j] + (dt/(dx*dx))*( v*(h[xx][j] + h[x][j] - 2*h[i][j]) + v*(h[i][yy] + h[i][y] - 2*h[i][j]) + 0.125*lambda*(pow(h[x][j] - h[i][j],2) + pow(h[i][j] - h[xx][j],2) + (h[x][j] - h[i][j])*(h[i][j] - h[xx][j])) + 0.125*lambda*(pow(h[i][y] - h[i][j],2) + pow(h[i][j] - h[i][yy],2) + (h[i][y] - h[i][j])*(h[i][j] - h[i][yy]))  ) + sigma*sqrt(12*dt)*RandomDouble(-0.5,0.5);
				
			}
		}

		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < L; ++j)
			{
				h[i][j] = h_temp[i][j];
				avg += h[i][j];
			}
			
		}

		avg = avg/(L*L);  //mean surface height

		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < L; ++j)
			{
				rms += pow(avg - h[i][j],2);
			}
		}
		rms = sqrt((rms/(L*L)));

		rough << t << "," << rms << endl;
	}
	rough.close();

	ofstream output("KPZ_2d_noise.csv");
	for (int i = 0; i < L; ++i)
	{
		for (int j = 0; j < L; ++j)
		{
			if (j != L-1)
			{
				output << h[i][j] << ",";
			} 
			else{
				output << h[i][j];
			}
		}
		output << endl;
	}
	output.close();
	return 0;
}