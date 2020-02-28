// Full 3d Monte Carlo Simulation ...

//** 

//** Chris Lehnen, Fourth Edition (5/28/10) **
//   Removed unused functions, streamlined code, and fixed bugs.
//   Implemented new diffusion methods: Modified Coordination method with Mullins statistics.
//                                      Escape method with Boltzmann statistics.
//   Implemented pause, resume, save, and load functionality.
//   Fully parameterized code - no direct modifications nessesary for simulaiton runs.
//   Added connical uniform flux distribution.
//   Separated flux distribution from substrate orientation.
//   Reworked re-emission to include specular reflection, cosine emission, and uniform hemispherical emission.
//   Added custom template option to load a template off of a text file.
//   Implemented dynamic memory control on lattice arrays.
//   Added .description file for use as parameter reference and a log file.


//** Ryan Badeau, Third Edition (1/16/09) **
//   Removed Outdated Code. Corrected Probability Distribution.
//   Added Flat Template and Comments
//   Implemented particle re-emissiom

//** Matthew Pelliccione, Second Edition (01/24/07) **
//   The 3D lattice saving and loading functions have the changes file = fopen (fileName, "wb");
//   and file = fopen (fileName, "rb");  instead of "w" and "r" only for binary data I/O

//** Dexian Ye, First Edition (03/24/04):
//   Atomic shadowing effect is included
//   -If the particle is about the head on an occupied site at the next step of movement
//   or passes by an occupied site ("atomic shadowing") at the orthogonal six nearest
//   directions (up,down, and four sides) then deposit at that current position

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "vector2.cpp"
#include "rng.cpp"
//#include <conio.h>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cin;
using std::endl;

const double pi = 3.14159265;
const double kB = 0.00008617343;
clock_t start;
ofstream fout;

/*/Code Overview
The code uses a periodic lattice. This means that regular boundary conditions are enforced - if a particle flies off a side, it is reintroduced with its exiting velocity
on the other side. The lattice is a great simplification for time and memory reasons - individual 'particles' are treated as single lattice elements, and are essentially unit squares.
IMPORTANT: The lattice is saved as a char bit field to reduce memory consumption, which means that the char data type is masked.
/*/
int counter = 0;
// **** ** User Defined Parameter Block ** ****
//****************************************************************************************************************
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      ***************SIMULATION SIZE***************

int N = 1024;      //** Volume is NxNxL ...// L must be an integer multiple of "8"
int L = 1024;      //** This is due to bit formatting, used for memory efficiency.

long int T = 10000;    //** T gives the lifetime of the simulation. Remember that the actual number of evolution steps is T*R.
long int R = 10000;    //** R defines the time scale of the simulation, by defining the number of particles/T

//      ***********PHYSICAL BEHAVIOR***********

float P[10] = {.01, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 1};   //** P[n] = Prob. of nth order particle sticking ...
int emitDist = 1;             //** Determines how re-emitted particles will interact with the surface
                              //**    1 = Reflection. The particle will simply reflect off of the surface in a specular manner.
                              //**    2 = Cosine distribution. Particles will leave the surface with a flux density ~ cos(theta), where theta
                              //**        is the angle off of the suface normal.
                              //**    3 = Uniform distribution. Particles will leave the surface with a flux density that is uniform over the
                              //**        hemisphere centered over the surface normal.

bool ETCH = false;            //** Parameter ETCH defines the nature of the incident particle. If deposition is selected (ETCH=false), a new incident particle
                              //**    will be created, move ballistically until it strikes the surface, where it will deposit following a side-wall-sticking behavior.
                              //**    Etching will similarly introduce a ballistic particle, but when it strikes the surface, will remove a particle.
bool SoS = true;              //** Sets the simulation to run in Solid-on-Solid mode. True = Solid-on-Solid. False = Ballistic Aggregation
                              //** Solid-on-Solid denotes behavior where a particle, upon deposition, will fall vertically to the lowest available location.
int DiffType = 1;             //** Determines the diffusion type that will be used
                              //**    1 = Boltzman Escape. Particles use Boltzmann statistics to see if they have enough energy to
                              //**        break all bonds at current position and move to a random neighboring position.
                              //**    2 = Linear Escape. Particles use a linear relationship with respect to current coordination to determining
                              //**        the likelihood of moving to a random neighboring position.
                              //**    3 = Sighted Coord. Particles look at the change in coordination caused by a given movement.
                              //**        This coordination change is then used to determine the likelihood of movement.
                              //**    4 = Sighted Height. Particles look at the change in height caused by a given movement, and use
                              //**        this change to determine the likelihood of movement. Note: In the code, the variables used are the same as those
                              //**        for coordination, and are named as such, but they are actually measureing height.
int side = 5;                 //** Range from incident particle where diffusion will take place. Diffusion box side = 2*side + 1
int skips = 0;                //** Number of diffusion steps per incident particle
int validLim = 2;             //** Minimum number of neighbors a particle can have and still be considered "stable"
float sideDepProb = .7;       //** The likelihood of  an incoming to deposit when there are particles at its sides.
float farCoord = 0;           //** Contribution to a location's coordination given by particles 2 spaces away

//**    Diffusion Parameters: Sighted Coordination (Mullins) & Sighted Height (Edwards-Wilkinson)
//              (parantheses indicate the anticipated statistical behavior created by this method)
float critCoord = 10;         //** Coordination differece where movement likelihood reaches 100%
float flatChance = .20;       //** Movement likelihood when the coordination difference = 0 (1=100%)

//**    Diffusion Parameters: Boltzmann Escape Coordination (Unknown. Mullins?)
float En = .05;               //** Bond energy for each neighboor (ie. unit coordination), in eV
float Ea = -3;                //** Static activation energy added to the neighboor bond energy, in units of coordination
int Ts = 300;                 //** Temperature of substrate, in Kelvin

//**    Diffusion Parameters: Linear Escape Coordiantion (Mullins)
float coord100 = 1;           //** Coordination at which a particle is guarenteed to diffuse
float coord0 = 17;            //** Coordination at which a particle will never diffuse

//      ***********DISTRIBUTION AND ROTATION**********

//**    Incoming Particle Distribution Parameters
int dist = 1;                                 //** Parameter dist identifies the type of incident particle flux to be simulated.
                                              //** 1 = uniform, 2 = cosine, 3 = reversed cosine 4 = cospower
float cospower = (2.0);                       //** 
                                              //** Notes: Cosine distribution is a good approximation to incident flux in CVD and Sputter conditions
float angleInner = 0;                         //** Starting angle of the spread of particles in degrees. Used only during uniform distribution.
float angleOuter = 0;                         //** Ending angle of the spread of particles in degrees. Used only during uniform distribution.

int additional_height = 0;                    //**Additional height above substrate for particle to start at
//**    Substrate Orientation Parameters
float OblqAngle = 0;                          //** OblqAngle is the oblique angle in degrees that the axis of the flux distribution makes with the surface normal

//Rotation
bool rotation = false;				          //** Rotate the substrate about its normal if true
unsigned long int PartPerRot = 500000;        //** Number of particles per rotation. (1/PartPerRot) is the frequency of the sample rotation.
int direct = 1;                               //** Rotation direction. 1 = clockwise, -1 = counterclockwise

//Rotation Swing
bool rotswing = false;                        //** Oscillate the substrate about its normal if true. Oscillations are done in a sinusoidal fashion
                                              //** RotSwing and Rotation are superpositional
float RotAngle = 90;						  //** Angle in degrees through which the substrate rotates
unsigned long int PartPerRotSwing = 200000;   //** Number of particle per rot swing cycle.
float phiphase = 0;                           //** Phase offset for rotswing in degrees

//Tilt Swing
bool tiltswing = false;			              //** Tilt the sample to the side in an oscillatory manner if true. Oscillations are done sinusoidally
                                              //** TiltSwing and OblqAngle are superpositional.
float TiltAngle = 20;                         //** Angle in degrees through which the substrate tilts.
unsigned long int PartPerTiltSwing = 250000;  //** Number of particles per tilt swing. (1/nPerSwing) is the frequency of the swinging.
float thetaphase = 0;                         //** Phase offset for tiltswing in degrees

//         *****************TEMPLATE PARAMETERS*********************

//The following parameter block is for allocating the size and morphology of a prexisting template for deposition.
//For example, one can form a uniformly distributed array of plugs, hemispheres, cones, etc...


int UseTemplate = 0;                    //Template to be used. 0= flat, 1= square, 2= flat conical, 3= colloid, 4= one, 5= colloid

// Form the flat top connical template (~tungsten plugs):
int rbottom = 16;				        //** Bottom radius of the cone
int hFCone = 20+8;			            //** Height of the cone (+8 is for most of the time the program starts from an initial height 8)
int rtop = 16;                          //** Top radius of the cone
int dPeakToPeak = 128+1;	            //** Peak-to-peak (ie. center to center) distance of the cones measured from the centers

//Form the colloid template: (also uses the above)
int Rcoll = 5;      		            //** Radius of the colloid

//Form the square template:
int Hsq = 32 + 8;			            //** Height of square template (+8 is for most of the time the program starts from an initial height 8)
int Lsq = 4;					        //** Width of square template
int Dsq = 30;					        //** Seperation between squares

//            **************OUTPUT FILE PARAMETERS*******************

char baseName [100] = "Name_Holder";      //** baseName gives the output filename to be preappended to all files.
bool sl = true;                           //** Saves the 3D lattice files if entered true
//**    Saves the following lattice view files if entered true
bool sly = true;            //** y-axis slice
bool slx = true;            //** x-axis slice
bool slz = true;            //** top view
bool spd = false;

// Isoflux Parameters
bool highp = false;
int spawn_height=20;
int spawn_std = 10;
int kill_height=200;
int killHDiff = 20;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//******************************************************************************************************
//  *** ** End of User Defined Parameters ** ***

//  *** ** Global Variables ** ***
//*********************************************************************************************************
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       *************SIMULATION STATE*****************

unsigned char ** Surface;      //** Main array containing the surface structure
particle Particle;             //** Particle object
int * H;                       //** Single valued height array
int ** Hcs;                    //** Single valued cross sectional arrays
int ** Hcs2;

unsigned long int nmov = 1;    //** Variable for moving particle, keeps track of how many steps the particle has moved from initialization

long int h0;                   //** Base height offset. Used when the Surface array data must be shifted up or down
long int Min;                  //** Minimum surface height as measured from a top-down perspective. Essentially min(H)
long int Max;                  //** Maximun surface height as measured from a top-down perspective. Essentially max(H)
long int MinNum;               //** Number of locations with height = Min
long int MaxNum;               //** Number of locations with height = Max

unsigned long int RotCount = 0;         //** Variable for keeping track of substrate rotation
unsigned long int RotSwingCount = 0;    //** Variable for keeping track of substrate rot swing
unsigned long int TiltSwingCount = 0;   //** Variable for keeping track of substrate tilt swing

//      *************WORKING STORAGE****************
int XTable[3*3*3];     //** Tables to keep track of a particle's neighbors for the recovery function
int YTable[3*3*3];
int ZTable[3*3*3];
double CDiff[3];       //** Tables for probability generation used in diffuse function
double CCoif[2];
double CConst[2];

// Table for fast and easy indexing through neighbor positions
const int S[26][3] = {{0,0,1},{0,0,-1},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},                           //** Nearest neighbors (0,5)
                    {1,0,-1},{-1,0,-1},{0,1,-1},{0,-1,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},      //** Next Nearest neighbors (6,17) ... SoS Nearest and Next Nearest neighbors (6,13)
                    {1,0,1},{-1,0,1},{0,1,1},{0,-1,1},
                    {1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1},{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1}};   //** Next Next Nearest neighbors (18,25)

//       ****************DATA RECORDING******************
int SaveTime[51] = {1, 2, 4, 6, 10,                               //** Determines when lattice images are taken. The values are such that when placed on a log scale
                    16, 25, 40, 63, 100,                          //** the images will appear to have been taken on roughly linear intervals
                    158, 251, 398, 631, 1000,
                    1250, 1500, 1750, 2000, 
                    2250, 2500, 2750, 3000, 
                    3250, 3500, 3750, 4000, 
                    4250, 4500, 4750, 5000, 
                    5250, 5500, 5750, 6000, 
                    6250, 6500, 6750, 7000, 
                    7250, 7500, 7750, 8000, 
                    8250, 8500, 8750, 9000, 
                    9250, 9500, 9750, 10000, };
		    
int ** Porosity;
                    
double HeightSave[10]={11, 13, 48, 75, 154, 500};
int SaveNumber = 0;
int HeightNumber=0;                                               //** Index number for SaveTime

// The following are used to record various statistical measures during the simulation
long int numCoord = 0;
double avgCoord = 0;
double preCoord = 0;
double rmsCoord = 0;
long int skipcount = 0;
double avgFinal = 0;
double preFinal = 0;
double rmsFinal = 0;
double avgDown = 0;
double avgUp = 0;
long int foundcount = 0;
long int recovered = 0;
long int deleted = 0;
double height = 0;
double width = 0;
double slope = 0;

bool WentOutside = false;           //** If turned to true during the simulation, it indicates that some part of the top surface
                                    //** of the film went outside of the Surface array at some point during the simulation.

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************************************************
//  *** ** End of Global Variables ** ***

//  *** ** Function Definitions** ***
//*********************************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//        *************SIM STATE ACCESS FUNCTIONS****************

// Gets the surface height at (i,j)
// Note: h(i,j) can be read or written
inline int & h (int i, int j)
{
   while (i < 0)
      i = i + N;
   while (i >= N)
      i = i - N;
   while (j < 0)
      j = j + N;
   while (j >= N)
      j = j - N;

   return H[j*N + i];
}

// Bit function to check if a site is occupied. Return true if site x,y,z is occupied,
// and false if unoccupied.
inline bool occupied (int x, int y, int z)
{
   while (x < 0)
      x = x + N;
   while (x >= N)
      x = x - N;
   while (y < 0)
      y = y + N;
   while (y >= N)
      y = y - N;

   if (z >= 0)
      return ((Surface[y*N + x][z/8] & (1 << (z%8))) != 0);
   else                         //** If the location is below the lattice, then only single
      return (H[y*N + x] >= z); //** valued data exists. So the position is considered
                                //** occupied if the column height is greater than the posisiton.
}

// Fills in the height array H
void findH ()
{
    for (int j = 0;j < N;++j)
   {
      for (int i = 0;i < N;++i)
      {
         H[j*N + i] = 0;

         for (int k = L - 1;k >= 0;--k)
         {
            if (occupied (i,j,k))
            {
               H[j*N + i] = k;
               break;
            }
         }
      }
   }
}

void findP () {
    int occu = 0;
    for (int k = 0; k < N; k++) {
        occu = 0; 
        for (int j = 0;j < N;++j)
        {
	    for (int i = 0;i < N;++i)
	    {
	        if (occupied (i,j,k))
	        {
		  occu++;
	        }
	    }
	}
	Porosity[0][k] = k;
	Porosity[1][k] = occu;
	Porosity[2][k] = N*N - occu;
    }
}

// Finds the cross-sectional depth data for Hcs
void findHcs ()
{
   int i = (int) (N/2);

   for (int j = 0;j < N;++j)
   {
      for (int k = 0;k < L;++k)
      {
         Hcs[j][k] = 0;

         for (int i = int(N/2);i >= 0;--i)
         {
            if (occupied (i,j,k))
            {
               Hcs[j][k] = i;

               break;
            }
         }
      }
   }
}

// Finds the cross-sectional depth data for Hcs2
void findHcs2 ()
{
   int j = (int) (N/2);

   for (int i = 0;i < N;++i)
   {
      for (int k = 0;k < L;++k)
      {
         Hcs2[i][k] = 0;

         for (int j = int(N/2);j >= 0;--j)
         {
            if (occupied (i,j,k))
            {
               Hcs2[i][k] = j;

               break;
            }
         }
      }
   }
}

//The findMin and findMax functions help save time by identifying surface boundaries,
//and are useful for functional definitions of surface (ie, when overhang is ignored).
void findMin ()
{
   Min = H[0];

   MinNum = 0;

   for (int j = 0;j < N;++j)
   {
      for (int i = 0;i < N;++i)
      {
         if (H[j*N + i] < Min)
         {
            Min = H[j*N + i];
 	        MinNum = 1;
         }
	     else if (H[j*N + i] == Min)
	     {
            ++MinNum;
         }
      }
   }
}

void findMax ()
{
   Max = H[0];

   MaxNum = 0;

   for (int j = 0;j < N;++j)
   {
      for (int i = 0;i < N;++i)
      {
         if (H[j*N + i] > Max)
         {
            Max = H[j*N + i];
	        MaxNum = 1;
         }
	     else if (H[j*N + i] == Max)
	     {
            ++MaxNum;
         }
      }
   }
}

// Function for finding the normal to a surface, regardless of orientation.
// It works under the principle that the final normal vector is the sum of the
// contributions from all neighboring particles.
inline vector2 normal(int x, int y, int z)
{
    vector2 vec = zeroVector;
    int i;

    // Look at all nearest neighbors, and add to the normal the vectors opposite
    // of the respective directions from the space in question.
    for (i = 0; i < 6; ++i)
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
        {
            vec.x += -S[i][0];
            vec.y += -S[i][1];
            vec.z += -S[i][2];
        }
    }
    // Do the same for the next nearest neighbors, but scale the contributions
    // so that the contributing vector still has a length of 1.
    for (i = 6; i < 18; ++i)
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
        {
            vec.x += -S[i][0]/sqrt(2);
            vec.y += -S[i][1]/sqrt(2);
            vec.z += -S[i][2]/sqrt(2);
        }
    }
    // Likewise for next next nearest.
    for (i = 18; i < 26; ++i)
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
        {
            vec.x += -S[i][0]/sqrt(3);
            vec.y += -S[i][1]/sqrt(3);
            vec.z += -S[i][2]/sqrt(3);
        }
    }
    if (norm(vec) == 0)
        return vec;

    return vec/norm(vec);   //** Scale the normal vector to be of unit length.
}

// Find the coordination at location (x, y, z)
inline float coordination (const int x, const int y, const int z)
{
    if (DiffType == 4)  //** In sighted height diffusion, the coordination is counted simply as the height
        return z;

    int near = 0;

    //Near neighboors
    for (int i = 0; i < 26; ++i) //just nearest neighbors
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
            ++near;
    }

    if (farCoord <= 0.001)
       return near;

    int X, Y, Z;
    int far = 0;

    //Far neighboors
    for (X = x-2; X <= x+2; X+=4)
    {
        for (Y = y-2; Y <= y+2; ++Y)
        {
            for (Z = z-2; Z <= z+2; ++Z)
            {
                if (occupied(X, Y, Z))
                {
                    ++far;
                }
            }
        }
    }
    for (X = x-1; X <= x+1; ++X)
    {
        for (Y = y-2; Y <= y+2; Y+=4)
        {
            for (Z = z-2; Z <= z+2; ++Z)
            {
                if (occupied(X, Y, Z))
                {
                    ++far;
                }
            }
        }
    }
    for (X = x-1; X <= x+1; ++X)
    {
        for (Y = y-1; Y <= y+1; ++Y)
        {
            for (Z = z-2; Z <= z+2; Z+=4)
            {
                if (occupied(X, Y, Z))
                {
                    ++far;
                }
            }
        }
    }

    return (((float)near) + ((float)far)*farCoord);
}

// Determine if the location is valid to support a particle. That is, if it has the
// required number of total neighbors and at least one nearest neighbor. Requiring validity
// is nessesary to prevent pairs (validLim > 1) or triplets (validLim > 2) of particles
// from rising off away from the surface during diffusion.
inline bool valid (const int x, const int y, const int z)
{
    int Num = 0;
    int i;
    //Nearest neighboors
    for (i = 0; i < 6; ++i)
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
            ++Num;
        if (Num >= validLim)
            return true;
    }

    if (Num == 0)
       return false;

    //Other neighboors
    for (i = 6; i < 26; ++i)
    {
        if (occupied(x + S[i][0], y + S[i][1], z + S[i][2]))
        {
            ++Num;
            if (Num >= validLim)
                return true;
        }
    }

    return false;
}

//      ***********SIMULATION ACTION FUNCTIONS******************
//

//Call function to etch (remove) a particle at location (i,j,k)
//****Note: This function has not been tested or reviewed since major program alterations
//    and is not guarenteed functional.
inline void etch (int i, int j, int k)
{
    int HH, kk;
                                 //** Use this function to modify the surface since it will modify
    while (i < 0)                //** the facilitative simulation state automatically ...
        i = i + N;
    while (i >= N)
        i = i - N;
    while (j < 0)
        j = j + N;
    while (j >= N)
        j = j - N;

    HH = H[j*N + i];

    if (k >= 0)
        Surface[j*N + i][k/8] = Surface[j*N + i][k/8] & (~(1 << (k%8)));


    if (k == HH)
    {
        for (kk = k;kk >= 0;--kk)
        {
            if (occupied (i,j,kk))
            {
                H[j*N + i] = kk;
                break;
            }
        }

        if (HH == Max)
        {
            --MaxNum;
            if (MaxNum == 0)
                findMax ();
        }

        if (kk < Min)
        {
            Min = kk;
            MinNum = 1;
            return;
        }

        if (kk == Min)
            ++MinNum;
    }
    return;
}

//Call function to deposit (add) a particle at location (i,j,k)
inline void deposit (int i, int j, int k)
{
   int HH;
                                 //** Use this function to modify the
   while (i < 0)                 //** surface since it will modify
      i = i + N;                 //** the facilitative simulation
   while (i >= N)                //** state automatically ...
      i = i - N;
   while (j < 0)
      j = j + N;
   while (j >= N)
      j = j - N;

   HH = H[j*N + i];

   if (k >= 0)
      Surface[j*N + i][k/8] = Surface[j*N + i][k/8] | (1 << (k%8));

   if (k > HH)
      H[j*N + i] = k;
   else
      return;

   if (HH == Min)
   {
      --MinNum;
	  if (MinNum == 0)
         findMin ();
   }

   if (k > Max)
   {
      Max = k;
      MaxNum = 1;
      return;
   }

   if (k == Max)
      ++MaxNum;

   return;
}
// Create a new particle at a random position above the surface maximum with a
// movement direction in accordance with the distribution settings.
// This creates a particle for low pressure deposition.
inline void initParticleLow (particle & P)
{
    P.xn = RandomInt(0,N-1);
    P.yn = RandomInt(0,N-1);
    P.zn = Max + 1 + additional_height;

    P.x = P.xn;
    P.y = P.yn;
    P.z = P.zn;

    double theta;
    double phi;

    P.vz = 0;
    while (P.vz > -0.01)   //** Make sure the particle has reasonable downward motion as to not stall program. This upper limit corresponds to an oblique angle of 89.4 degrees.
    {
        if(dist == 1) //** Uniform spread flux: Particles are evenly distributed throughout a cone (or hollow cone) with angleInner <= theta <= angleOuter
        {
            theta = acos(RandomDouble(cos(angleOuter*pi/180), cos(angleInner*pi/180)));
            phi = RandomDouble(0, 2*pi);
        }
        else //** Cosine distribution: Particles are distributed with density ~ cos(theta)
        {
             //VERY IMPORTANT, to generate a correct cosine distributed theta, the function should be the arcsin(sqrt(x)), where x the random number.
             //Leaving out the sqrt will actually generate a cotangent distribution!
             //See the article by Greenwood, 2002, which carries out the polar derivation
             //The correct and incorrect generation of a cosine distribution of scattered particles for Monte-Carlo modelling of vacuum systems
            theta = asin(sqrt(RandomDouble(0,1)));
            phi = RandomDouble(0, 2*pi);
        }

        // Generate a unit directional vector based on the angular direction
        P.vz = (cos(theta));
        P.vx = (sin(theta)*cos(phi));
        P.vy = (sin(theta)*sin(phi));

        // Rotate directional vector to represent substrate orientation.
        //
        theta = OblqAngle*pi/180; //** Tilt the substrate by OblqAngle
        if(tiltswing) //** Add tiltswing to OblqAngle if activated
        {
            theta = theta + TiltAngle*(pi/180)*sin(float(TiltSwingCount)/PartPerTiltSwing*2*pi + thetaphase*(pi/180));
            ++TiltSwingCount;
            if (TiltSwingCount >= PartPerTiltSwing)
                TiltSwingCount -= PartPerTiltSwing;
        }
        if(rotation) //** Rotate the substrate about normal if activated
        {
            phi = float(RotCount)/PartPerRot*2*pi*direct;
            ++RotCount;
            if (RotCount >= PartPerRot)
                RotCount -= PartPerRot;
        }
        else  //** Otherwise reset phi
            phi = 0;
        if(rotswing) //** Add rotswing to rotation if activated
        {
            phi = phi + RotAngle*(pi/180)* sin(float(RotSwingCount)/PartPerRotSwing*2*pi + phiphase*(pi/180));
            ++RotSwingCount;
            if (RotSwingCount >= PartPerRotSwing)
                RotSwingCount -= PartPerRotSwing;
        }

        // Apply two rotation matrices to get the final directional vector
        // First the vector is rotated around the y-axis by theta, then the result is rotated around the z-axis by phi
        double tempx = P.vx*cos(phi)*cos(theta) - P.vy*sin(phi) + P.vz*cos(phi)*sin(theta);
        double tempy = P.vx*sin(phi)*cos(theta) + P.vy*cos(phi) + P.vz*sin(phi)*sin(theta);
        double tempz = -P.vx*sin(theta) + P.vz*cos(theta);
        P.vx = tempx;
        P.vy = tempy;
        P.vz = -tempz; //** The particle is made to go downwards here
    }
    //cout << P.x << " " << P.y << " " << P.z << " " << theta << " " << phi << " " << P.vx << " " << P.vy << " " << P.vz << endl; // To check initialization position, angles, and direction
    P.impact = 0;
    setPInternal (P);
    nmov = 1;
    return;
}

// Create a new particle at a random position above the surface 
// The velocities are distributed randomly over a sphere
// This creates a particle for High pressure deposition. (or simply conformal) 
inline void initParticleHighPressure (particle & P) {
// 
    double theta;
    int tempi = RandomInt(0,N-1);
    int tempj = RandomInt(0,N-1);
    P.xn = tempi;
    P.yn = tempj;
    int localHeight = h(tempi,tempj);
    int added_height = ceil(RandomGaussian(((double)spawn_height), ((double)spawn_std)));
    if (added_height < 2) 
    {
        added_height = 2;
    }
    
    P.zn = localHeight + added_height; // Changed to simple random column choice with some hight above the current column the one added can be changed to anything.
    kill_height = localHeight + killHDiff;
    P.x = P.xn;
    P.y = P.yn;
    P.z = P.zn;
    if (dist == 4) {        
        theta = acos(pow((1 - RandomDouble(0,1)), 1.0/(cospower + 1)));        
    } 
    else if (dist == 2) {
        
        theta = asin(sqrt(RandomDouble(0,1)));
    }
    else 
    {
        
        theta = acos(RandomDouble(-1,1));
    } 
    
    double phi = RandomDouble(0, 2*pi);

    P.vz = (cos(theta));
    P.vx = (sin(theta)*cos(phi));
    P.vy = (sin(theta)*sin(phi));
    
    if (dist != 1) {
        P.vz = -P.vz;
        
    }

    //cout << P.x << " " << P.y << " " << P.z << " " << theta << " " << phi << " " << P.vx << " " << P.vy << " " << P.vz << endl; // To check initialization position, angles, and direction
    P.impact = 0;
    setPInternal (P);
    nmov = 1;
    return;
}



// Function for particle reemission. It finds the normal of the surface at the particle's deposition
// location and gives the particle a new direction depending on this normal and the selected options.
inline void emitParticle (particle & P)
{
    //printf("a");
    vector2 normvec = normal(P.xl, P.yl, P.zl);   //** Get the normal vector.
    vector2 partvec = {-P.vx, -P.vy, -P.vz};      //** Get the vector on the same line as the particle's direction, but away from the surface.

    if (emitDist == 1)                              //** Specular Reflection. New direction is a reflection off of the surface.
    {
        //printf("b");
        partvec = reflected(normvec, partvec);
    }
    else
    {
        double theta;
        double phi = RandomDouble(0, 2*pi);
        if (emitDist == 2)                          //** Cosine Spread. Like the deposition distribution, but centered along the normal.
            theta = asin(sqrt(RandomDouble(0,1)));
        else                                        //** Uniform Spread . A uniform distribution over a hemisphere centered along the normal.
            theta = acos(RandomDouble(0,1));

        // Initial distribution
        P.vz = (cos(theta));
        P.vx = (sin(theta)*cos(phi));
        P.vy = (sin(theta)*sin(phi));

        // Find the needed rotation angles to line up the distribution with the normal.
        theta = acos(normvec.z);
        if (normvec.x == 0)
        {
            if (normvec.y > 0)
                phi = pi/2;
            else
                phi = 3*pi/2;
        }
        else
        {
            phi = atan(normvec.y/normvec.x);
            if (normvec.x < 0)
               phi += pi;
        }

        // Apply two rotation matrices to get the final directional vector
        // First the vector is rotated around the y-axis by theta, then the result is rotated around the z-axis by phi
        // This should align the axis of the distribution with the normal of the surface.
        partvec.x = P.vx*cos(phi)*cos(theta) - P.vy*sin(phi) + P.vz*cos(phi)*sin(theta);
        partvec.y = P.vx*sin(phi)*cos(theta) + P.vy*cos(phi) + P.vz*sin(phi)*sin(theta);
        partvec.z = -P.vx*sin(theta) + P.vz*cos(theta);
    }

    //printf("c");
    P.vx = partvec.x;
    P.vy = partvec.y;
    P.vz = partvec.z;

    P.x = P.xl;
    P.y = P.yl;
    P.z = P.zl;

    P.xn = P.x;
    P.yn = P.y;
    P.zn = P.z;

    setPInternal (P);
    nmov = 1;
    ++P.impact;
    //printf("d ");
}

// Moves the particle one step in the appropriate direction, enforcing periodic boundary conditions.
inline void moveParticle (particle & P)
{
   // Make the current position be the "last" position
   P.xl = P.x;
   P.yl = P.y;
   P.zl = P.z;

   // Update the current position to be one step farther
   P.x = (int)(round(P.xn + P.vx*nmov));
   P.y = (int)(round(P.yn + P.vy*nmov));

   // Enforce boundary conditions
   while(P.x < 0)
       P.x += N;
   while(P.x >= N)
       P.x -= N;
   while(P.y < 0)
       P.y += N;
   while(P.y >= N)
       P.y -= N;

   P.z = (int)(round(P.zn + P.vz*nmov));

   ++nmov; //** Increase step count for next move
   return;
}

// This function takes in the final coordination value (though depending on the method it may
// actually be a coordination difference or height difference) and determines if the particle will
// move or not by comparing a calculated probability with a random number.
inline bool diffuse(float C)
{
    if (DiffType == 1)
    {
        return (RandomDouble(0,1) < exp(-(En*C+Ea)/(kB*Ts)) ); //** Boltzman distribution
    }
    if (DiffType == 2)
    {
        return (RandomDouble(0,1) < (coord0-C)/(coord0-coord100) ); //** Linear probability
    }
    if (DiffType == 3 || DiffType == 4) //** Piecewise linear function that may allow for non-zero probability despite
    {                                   //** negative coordination differences (rises in energy). However the
        if (C < CDiff[0])               //** function is such that the net movement probability in the upward
            return false;               //** direction is always directly proportional to the rise in coordination.
        if (C < CDiff[1])
            return (RandomDouble(0,1) < (CCoif[0]*C + CConst[0]));
        if (C < CDiff[2])
            return (RandomDouble(0,1) < (CCoif[1]*C + CConst[1]));
        return true;
    }
    return false;
}

inline void checklocal (int x0, int y0, int z0);

// This function deals with the situation where a particle is deposited at an unstable location
// (which can happen in instances of direct impact) or a particle is residing in a location that
// becomes unstable due to the removal of a neighbor. Upon such an event, the particle must be moved
// to a neighboring, stable, location.
inline void recover(int x, int y, int z, bool check)
{
    int n, i;
    int mn = 0;

    ++foundcount;

    if (SoS)            //** In solid-on-solid mode, only the top particles are of concern
        z = h(x, y);

    etch(x, y, z);      //** Particle is temporariliy removed so as not to interfere with coordination measurements

    // Find all valid neighboring locations
    if (SoS)
    {
        for (i = 6; i < 14; ++i)
        {
            if (valid(x + S[i][0], y + S[i][1], h(x+S[i][0], y+S[i][1]) + 1 ))
            {
                XTable[mn] = x + S[i][0];
                YTable[mn] = y + S[i][1];
                ZTable[mn] = h(x+S[i][0], y+S[i][1]) + 1;
                ++mn;
            }
        }
    }
    else
    {
        for (i = 0; i < 26; ++i)
        {
            if ((!occupied(x + S[i][0], y + S[i][1], z + S[i][2])) && valid(x + S[i][0], y + S[i][1], z + S[i][2]))
            {
                XTable[mn] = x + S[i][0];
                YTable[mn] = y + S[i][1];
                ZTable[mn] = z + S[i][2];
                ++mn;
            }
        }
    }
    // If a sutable location cannot be found, the particle is assumed to
    // fly off the surface and is deleted (or just not put back).
    if (mn == 0)
    {
       ++deleted;
       if (check)
           checklocal(x, y, z);  //** If a particle is removed then the surrounding particles must be checked to ensure
                                                //** they are still valid, except in the case of a newly deposited particle.
       return;
    }

    ++recovered;
    // If possible, the process of choosing a final movement location is done in a manner
    // consistant with the current diffusion method.
    if (DiffType == 1 || DiffType == 2)
    {
        // In these two methods, the particles have no knowledge about the conditions of their destination,
        // so there is no push towards one location versus another
        n = RandomInt(0, mn-1);
        deposit(XTable[n], YTable[n], ZTable[n]);
        if (check)
            checklocal(x, y, z);
    }
    else if (DiffType == 3 || DiffType == 4)
    {
        // For these two methods, a random location from the list of valid locations is chosen, then the
        // standard diffusion procedure is done for that location. If the particle moves, the process is
        // finished. If not, a new random location is chosen. This continues until the particle moves
        // or there is only one possible location left, at which point the particle moves to that location.
        float cc = coordination(x, y, z);
        while (mn > 1)
        {
            n = RandomInt(0, mn-1);
            if (diffuse(coordination(XTable[n], YTable[n], ZTable[n]) - cc))
            {
                deposit(XTable[n], YTable[n], ZTable[n]);
                if (check)
                    checklocal(x, y, z);
                return;
            }
            else
            {
                XTable[n] = XTable[mn-1];
                YTable[n] = YTable[mn-1];
                ZTable[n] = ZTable[mn-1];
                --mn;
            }
        }

        deposit(XTable[0], YTable[0], ZTable[0]);
        if (check)
            checklocal(x, y, z);
    }

    return;
}

// This function checks all neighbors of a given location and makes sure that they are stable.
// If a particle is present at an unstable location, then the recover function is called on that
// location. The recover function will then call this function again around the particle that
// moved. Hence the whole process is recursive, and will cascade though a series of instabilities
// until all particles have obtained stable locations.
inline void checklocal (int x0, int y0, int z0)
{
    int xc, yc, zc;
    for (int i = 0; i < 26; ++i)
    {
        xc = x0 + S[i][0];
        yc = y0 + S[i][1];
        zc = z0 + S[i][2];
        if (occupied(xc, yc, zc) && (!valid(xc, yc, zc)))
            recover(xc, yc, zc, true);
    }

    return;
}

// This function controls the diffusion behavior of a particle. The basic process is as follows:
// 1) A random direction is chosen. 2) If the location in that direction is valid and unoccupied
// then coordination values are measured depending on the selected diffusion type. 3) The final
// coordination value is used to decide if the particle will move, and if so, the particle is
// moved 4) If any part of the above process failed, the function is exited, and the simulation
// continues.
// The purpose of only attempting one direction per diffusion step is to ensure a uniform
// distribution of movement attempts so that the actual movement behavior is dependent only on the
// statistics of the selected diffusion type and not on the relative number/position of
// neighboring particles.
#include "diffusion.h"
inline void skip (int X1, int Y1, int Z1)
{
    int n, X2, Y2, Z2;
    unsigned char mask;
    int origH;

    if (SoS)                      //** In solid-on-solid mode, only top particles are considered for diffusion
    {
        Z1 = h(X1, Y1);
        n = RandomInt(6, 13);
        X2 = X1 + S[n][0];
        Y2 = Y1 + S[n][1];
        Z2 = h(X2, Y2) + 1;
    }
    else
    {
        if (!occupied(X1, Y1, Z1))   //** Make sure that there is actually a particle to move
            return;

        n = RandomInt(0, 25);        //** Choose a random direction and final movement location
        X2 = X1 + S[n][0];
        Y2 = Y1 + S[n][1];
        Z2 = Z1 + S[n][2];

        if (occupied(X2, Y2, Z2))    //** Make sure that the final location is unoccupied
            return;
    }

    // This whole block is used to check if the final location is valid.
    // It is first nessesary to remove the current particle so that it does not affect
    // the validity and coordination measurements. Basically, it makes sure the particle
    // does not view itself as a neighbor.
    while (X1 < 0)  X1 = X1 + N;
    while (Y1 < 0)  Y1 = Y1 + N;
    while (X1 >= N) X1 = X1 - N;
    while (Y1 >= N) Y1 = Y1 - N;

    if (Z1 >= 0)
    {
        mask = 1 << (Z1%8);
        Surface[Y1*N + X1][Z1/8] = Surface[Y1*N + X1][Z1/8] ^ mask;
    }
    else
    {
        origH = H[Y1*N + X1];    //** If the particle is below the lattice, adjust the height array instead.
        H[Y1*N + X1] = Z1 - 1;
    }

    if (!valid(X2, Y2, Z2))      //** Final location is checked for validity
    {
        if (Z1 >= 0)             //** If the validity check fails, return the particle to its position.
            Surface[Y1*N + X1][Z1/8] = Surface[Y1*N + X1][Z1/8] ^ mask;
        else
            H[Y1*N + X1] = origH;

        return;                  //** and leave the function now.
    }

    float coord;
    // Find the corresponding coordination value (or difference value) depending on diffusion type
    if (DiffType == 1 || DiffType == 2)
        coord = coordination(X1, Y1, Z1);
    else if (DiffType == 3)
        coord = coordination(X2, Y2, Z2) - coordination(X1, Y1, Z1);
    else if (DiffType == 4)
        coord = coordination(X1, Y1, Z1) - coordination(X2, Y2, Z2);

    if (Z1 >= 0)
        Surface[Y1*N + X1][Z1/8] = Surface[Y1*N + X1][Z1/8] ^ mask; //** Return particle to its position
    else
        H[Y1*N + X1] = origH;

    // Statistics collection
    avgCoord += coord;
    rmsCoord += (coord - preCoord)*(coord - preCoord);
    ++numCoord;

    if (diffuse(coord)) //** Decide if the particle will move
    {
        etch(X1, Y1, Z1);           //** Move the particle. It is nessesary to use the etch function
        deposit(X2, Y2, Z2);        //** as corresponing height data will be updated.

        // Statistics collection
        ++skipcount;
        avgFinal += coord;
        rmsFinal += (coord - preFinal)*(coord - preFinal);
        if (coord >= 0) ++avgDown;
        else ++avgUp;

        checklocal(X1, Y1, Z1);     //** Since a particle left this location, check to ensure the
                                    //** stability of neighbors.
    }
    return;
}

// Shifts all data in the lattice down by 32 units to make room for further growth
void ShiftDown ()
{
   for (int i = 0;i < N;++i)
   {
      for (int j = 0;j < N;++j)
      {
         H [i*N + j] = H [i*N + j] - 32;

         for (int k = 0;k < (L-32)/8;++k)
         {
            Surface [i*N + j][k] = Surface [i*N + j][k + 4];
         }
         for (int k = (L-32)/8;k < L/8;++k)
         {
            Surface [i*N + j][k] = 0;
         }
      }
   }

   Max = Max - 32;
   Min = Min - 32;
   if (Min < 0)            //** Checks to see if the surface of the film dipped below the lattice. Those regions of the surface
       WentOutside = true; //** are now treated as single-valued. Such an occurance could lead to erroneous behavior in the future
   h0 = h0 + 32;           //** of the simulation, so the fact that it occured is noted.
   return;
}

// Shifts all data in the lattice up by 32 units to make room for further etching
void ShiftUp ()
{
   for (int i = 0;i < N;++i)
   {
      for (int j = 0;j < N;++j)
      {
         H [i*N + j] = H [i*N + j] + 32;

         for (int k = L/8 - 1;k >= 4;--k)
         {
            Surface [i*N + j][k] = Surface [i*N + j][k - 4];
         }

         for ( int k = 0;k < 4;++k)
         {
            Surface [i*N + j][k] = 255;
         }
      }
   }

   Max = Max + 32;
   Min = Min + 32;
   if (Max > L-1)          //** Checks to see if the surface of the film rose above the lattice. Those regions of the surface
       WentOutside = true; //** have now been effectively cut off. Such an occurance could lead to erroneous behavior in the future
   h0 = h0 - 32;           //** of the simulation, so the fact that it occured is noted.
   return;
}

//     ***********SAVE DATA FUNCTIONS************

void saveSim (int n)
{
   int x, y, z;
   unsigned char value;

   char fileName[100];
   sprintf (fileName, "%s_Save.txt", baseName);
   FILE * file;
   file = fopen (fileName, "w");

    fprintf (file, "%i\n%i\n%il\n%i\n%il\n%i\n%il\n%fL\n%fL\n", (int) N, (int) L, (int) R, (int) n, (int) h0, (int) SaveNumber, (int) deleted, preCoord, preFinal);
    for (z=0; z<L/8; ++z)
    {
        for (y=0; y<N; ++y)
        {
            for (x=0; x<N; ++x)
            {
                value = Surface[y*N + x][z];
                fprintf (file, "%c", value);
            }
            fprintf (file, "\n");
        }
        fprintf (file, "\n");
   }


   fclose (file);
}

// Save the full 3d lattice
void saveLat (char * fileName)
{
   int x, y, z;
   unsigned char value;

   FILE * file;
   file = fopen (fileName, "w");

    for (z=0; z<L/8; ++z)
    {
        for (y=0; y<N; ++y)
        {
            for (x=0; x<N; ++x)
            {
                value = Surface[y*N + x][z];
                fprintf (file, "%c", value);
            }
            fprintf (file, "\n");
        }
        fprintf (file, "\n");
   }


   fclose (file);
}

// Save the H array as a text file for latter viewing as an image.
void saveHeight (char * fileName)
{
    FILE * file;
    file = fopen (fileName, "w");

    for (int i = 0;i < N;++i)
    {
        for (int j = 0;j < N;++j)
        {
            fprintf (file, "%i", (int) h(j,i));
            if (j < N - 1)
                fprintf (file, " ");
        }
        fprintf (file, "\n");
    }
    fclose (file);
}

void savePorosity (char * fileName) 
{
    FILE * file;
    file = fopen (fileName, "w");
    findP ();
    for (int i = 0;i < L;++i)
    {
        for (int j = 0;j < 3;++j)
        {
            fprintf (file, "%i", Porosity[j][i]);
            if (j < 2)
                fprintf (file, " ");
        }
        fprintf (file, "\n");
    }
    fclose (file);
    
}

// Save the Hcs array as a text file for latter viewing as an image.
void saveCrossSec3D (char * fileName)
{
    FILE * file;
    file = fopen (fileName, "w");

    for (int k = 0;k < L;++k)
    {
        for (int j = 0;j < N;++j)
        {
            fprintf (file, "%i", (int) Hcs[j][k]);
            if (j < N - 1) fprintf (file, " ");
        }
        fprintf (file, "\n");
    }

    // The viewing program stmCalc does not properly import non-square images.
    // For that reason, we must extend the image file above the top of the
    // lattice if L < N.
    for (int k = L;k < N;++k)
    {
        for (int j = 0;j < N;++j)
        {
            fprintf (file, "%i", 0);
            if (j < N - 1) fprintf (file, " ");
        }
        fprintf (file, "\n");
    }

    fclose (file);
}

// Save the Hcs2 array as a text file for latter viewing as an image.
void saveCrossSec3D2 (char * fileName)
{
    FILE * file;
    file = fopen (fileName, "w");

    for (int k = 0;k < L;++k)
    {
        for (int j = 0;j < N;++j)
        {
            fprintf (file, "%i", (int) Hcs2[j][k]);
            if (j < N - 1) fprintf (file, " ");
        }
        fprintf (file, "\n");
    }

    // See above.
    for (int k = L;k < N;++k)
    {
        for (int j = 0;j < N;++j)
        {
            fprintf (file, "%i", 0);
            if (j < N - 1) fprintf (file, " ");
        }
        fprintf (file, "\n");
    }

    fclose (file);
}

//        *********TEMPLATE FUNCTIONS************
//****Note (by Chris): The following template functions have not been checked or verified for lack of errors
//    for a long time. Neither have I studied them thoroughly enough to comment on their function. I would
//    advise reviewing their code first before use.
void formSquareTemplate ()
{
   //**Form the inital square template surface**
   int c=0;
   int label=1;
   int p =0 ;
   int pattern = 1;

   for (int i = 0; i<N; ++i) // check paterning in direction-i
   {
      ++p;
	  c = 0;
	  label=1;
	  if (pattern>0)
	  {
         for (int j = 0; j<N; ++j) //check patterning in direction-j
		 {
            ++c;
			if  (label>0)
			{
			   for (int k = 0; k<Hsq; ++k) { deposit(i,j,k); }
			   if (c==Lsq)
			   {
			      label = -1; //change patterning critaria
			      c = 0;
               }
			}
			if ( (label < 0) && (c==Dsq) )
			{
			   label = 1; //change patterning critaria
			   c = 0;
			}
         }

		 if (p==Lsq)
		 {
            pattern = -1; //change patterning critaria
			p = 0;
		 }
      }

      if ( (pattern<0) && (p==Dsq) )
	  {
         pattern = 1; //change patterning critaria
	  	 p = 0;
	  }
   }

   return;
}

void formColloidTemplate ()
{
   float hpow2;
   float heightcoll; //surface height of the spherical colloid
   int shiftRow = -1;
   int deltai;

   for (int i = 0; i<N; ++i)
   {
      shiftRow = -1;
      deltai = -(Rcoll/2);
      for (int j = 0; j<N; ++j)
      {
         //Spheres in bcc (100) order:
         //hpow2 = pow(Rcoll,2) - pow( (i%(2*Rcoll) -Rcoll),2) - pow( (j%(2*Rcoll) -Rcoll),2);

         //Spheres in closed pack (fcc (111)?) order:
         if((j%(2*Rcoll))==0)
         {
            shiftRow = shiftRow*(-1);
         }

         deltai = (Rcoll/2)*shiftRow;
         hpow2 = pow((double) Rcoll,2) - pow((double) ((i+deltai)%(2*Rcoll) -Rcoll),2) - pow((double) (j%(2*Rcoll) -Rcoll),2);

         if (hpow2<0) {hpow2 = 0;}
            heightcoll = sqrt(hpow2);

         for (int k = 0; k<(heightcoll+7); ++k)
         {
            deposit(i,j,k);
         }
      }
   }

   return;
}

void formFlatConicalTemplate ()
{
   int R = int((dPeakToPeak-1)/2);
   int r = rbottom;

   for(int k = 8; k<=hFCone; ++k) //k start at 8 since the substrate surface is at height=7.
   {
      r = int ( rtop + (rbottom-rtop)*((hFCone-7)-(k-7))/(hFCone-7) ); //the cone wall angle (angle between the wall and substrate plane) is invTan (hCone/(rbottom-rtop))

      for (int i = 0; i<N; ++i)
      {
         for (int j = 0; j<N; ++j)
         {
            //if  ( ( pow(( (i%2*R)-R),2) + pow(( (j%2*R)-R),2) ) <= pow(r,2) )
            //if  ( ( pow(( int(i/(2*R))-i),2) + pow(( int(j/(2*R))-j),2) ) <= pow(r,2) )
            if  ( ( pow((double)( (int((i+R)/(2*R)))*(2*R)-i),2) + pow((double)( (int((j+R)/(2*R)))*(2*R)-j),2) ) <= pow((double)r,2) )
            {
               deposit(i,j,k);
            }
         }
      }
   }

   return;
}

void formOne()
{
   int R = int((dPeakToPeak-1)/2);
   int r = rbottom;

   for(int k = hFCone; k<=hFCone; ++k) //k start at 8 since the substrate surface is at height=7.
   {
      //r = int ( rtop + (rbottom-rtop)*((hFCone-7)-(k-7))/(hFCone-7) ); //the cone wall angle (angle between the wall and substrate plane) is invTan (hCone/(rbottom-rtop))

      for (int i = N/2-r-1; i<N/2+r+1; ++i)
      {
         for (int j = N/2-r-1; j<N/2+r+1; ++j)
         {
            //if  ( ( pow(( (i%2*R)-R),2) + pow(( (j%2*R)-R),2) ) <= pow(r,2) )
            //if  ( ( pow(( int(i/(2*R))-i),2) + pow(( int(j/(2*R))-j),2) ) <= pow(r,2) )
            if  ( ( pow((double) (N/2-i),2) + pow((double)( N/2-j),2) ) <= pow((double)r,2) )
            {
			   deposit(i,j,k);
            }
         }
      }
   }

   return;
}

void formArcConicalTemplate ()
{
	int R = int((dPeakToPeak-1)/2);
	int r = rbottom;
    float rs;
	for(int k = 8; k<=hFCone+(rtop/2); ++k) //k start at 8 since the substrate surface is at height=7.
	{
	    if(k<=hFCone)
		{
            r = int ( rtop + (rbottom-rtop)*((hFCone-7)-(k-7))/(hFCone-7) ); //the cone wall angle
	        rs = pow((double)r,2);                          //(angle between the wall and substrate plane) is invTan (hCone/(rbottom-rtop))
        }
		else
		    rs = pow((double)(1.25*rtop),2) - pow((double)(0.75*rtop + k - hFCone),2);

		for (int i = 0; i<N; ++i)
		{
			for (int j = 0; j<N; ++j)
			{
				//if  ( ( pow(( (i%2*R)-R),2) + pow(( (j%2*R)-R),2) ) <= pow(r,2) )
				//if  ( ( pow(( int(i/(2*R))-i),2) + pow(( int(j/(2*R))-j),2) ) <= pow(r,2) )
				if  ( ( pow( (double)( (int((i+R)/(2*R)))*(2*R)-i),2) + pow( (double)( (int((j+R)/(2*R)))*(2*R)-j),2) ) <= rs )
				{
					deposit(i,j,k);
				}
			}
		}
    }
    return;
}

void formCustomTemplate()
{
    ifstream infile("template.z");
    if (!infile.good())                     //** Make sure the template file loaded properly.
    {
        cout << "Error: Template file does not exist.\nPress any key to close.\n";
        cin.get();
        infile.close();
        exit(1);
    }
    int value;
    for (int j=0; j<N; ++j)
    {
        for (int i=0; i<N; ++i)
        {
            infile >> value;
            if (!infile.good())             //** Make sure the end of the file has not been reached.
            {
                cout << "Error: End of template file reached.\nPress any key to close.\n";
                cin.get();
                infile.close();
                exit(1);
            }
            for (int k=0; k<value; ++k)     //** Build the surface up to the level of the template.
            {
                deposit(i,j,k);           //** 0 height on the template corresponds to 8 surface height.
            }
        }
    }
    infile >> value;
    if (infile.good())                      //** Make sure the end of the file has been reached, indicating the correct template size.
    {
        cout << "Error: Template larger than substrate. Template not read correctly.\nPress any key to close.\n";
        cin.get();
        infile.close();
        exit(1);
    }
    infile.close();
}


//        **************INITALIZATION FUNCTIONS****************

// Create a solid layer of particles 8 deep at the bottom of the lattice.
// Used as the initial surface for deposition simulations.
void initLow ()
{
   for (int i = 0;i < N;++i)
   {
      for (int j = 0;j < N;++j)
      {
         for (int k = 1;k < L/8;++k)
	     {
            Surface [j*N + i][k] = 0;
         }

         Surface [j*N + i][0] = 255;
      }
   }
}

// Fill the entire lattice with particles except for a layer at the top of the
// lattice 8 deep. Used as the initial surface for etching simulations.
void initHigh ()
{
   for (int i = 0;i < N;++i)
   {
      for (int j = 0;j < N;++j)
      {
         for (int k = 0;k < L/8 - 1;++k)
	     {
            Surface [j*N + i][k] = 255;
         }

         Surface [j*N + i][L/8 - 1] = 0;
      }
   }
}


// Initialization function run at the begining of every simulation. It allocates
// memory, creates the initial surface, etc.
int Initialize ()
{
    // Memory allocation.
    printf ("allocating memory ... \n");
    Porosity = new int * [3];
    Porosity[0] = new int [L];
    Porosity[1] = new int [L];
    Porosity[2] = new int [L];
    Surface = new unsigned char * [N*N];
    assert (Surface != 0);
    for (int i = 0;i < N*N;++i)
    {
        Surface [i] = new unsigned char [L/8];
        assert (Surface [i] != 0);
    }

    H = new int [N*N];
    Hcs = new int * [N];
    Hcs2 = new int * [N];
    for (int i = 0; i < N; ++i)
    {
        Hcs [i] = new int [L];
        Hcs2 [i] = new int [L];
    }

    // Check for save file
    printf ("checking for saved simulation\n");
    int loadn = 0;//LoadSave();

    if (loadn == 0)
    {
        // Initial surface creation
        printf ("initializing new surface and template...\n");

        if (ETCH == 1)
           initHigh ();
        else
           initLow ();

        // The flat template (stated as UseTemplate = 0) is really just the absence of a template,
        // so any value not from 1-5 will result in a flat surface.
        if (UseTemplate == 1)
            formSquareTemplate ();      //** Start with a templated surface formed by square-plugs
        else if (UseTemplate == 2)
            formFlatConicalTemplate();  //** Start with a templated surface formed by flat top cones.
        else if (UseTemplate == 3)
            formArcConicalTemplate();
        else if (UseTemplate == 4)
            formOne();
        else if (UseTemplate == 5)
            formColloidTemplate ();     //** Start with a templated surface formed by colloids
        else if (UseTemplate == 6)      //** Start with a surface defined by a custom template in the file "template.z"
            formCustomTemplate();
    }
    // Find initial surface data
    printf ("finding surface profile...\n");

    findH   ();
    findMin ();
    findMax ();
    if (loadn == 0)
        h0 = 0;

    // Initialize various other variables used in the simulation.
    printf ("initializing working variables...\n");
    WentOutside = false;
    if (loadn == 0)
    {
        deleted = 0;
        SaveNumber = 0;
    }

    // Initialize the tables that define the diffusion behavor for diffusion types 3 and 4.
    // The behavior is a series of linear relationships that form a piecewise curve extending from - to + critCoord and from 0 to 1 probability.
    // CDiff defines the bounds of the different regions
    // CCoif is the coeficient and CConst the constant in the equation P(x) = Ax + B
    if (critCoord <= 0.0001)    //** On/Off diffusion behavior where a particle is guarenteed to move
    {                           //** with coordination difference > 0 and prohibited when < 0
        CDiff[0] = -0.0001;
        CDiff[1] = 0.0001; CCoif[0] = 0; CConst[0] = flatChance;
        CDiff[2] = 0.0001; CCoif[1] = 0; CConst[1] = 2;
    }
    else
    {
        if (-critCoord <= -flatChance*2*critCoord)      //** Smoother diffusion behavior where down-coordination movements are allowed.
        {                                               //** Roughly, flatChance defines the shape of the behavior while critCoord defines the extent.
            CDiff[0] = -flatChance*2*critCoord;
            CDiff[1] = flatChance*2*critCoord; CCoif[0] = 1/(2*critCoord); CConst[0] = flatChance;
            CDiff[2] = critCoord; CCoif[1] = 1/critCoord; CConst[1] = 0;
        }
        else
        {
            CDiff[0] = -critCoord;
            CDiff[1] = (flatChance - 1)*2*critCoord; CCoif[0] = 1/critCoord; CConst[0] = 1;
            CDiff[2] = (1 - flatChance)*2*critCoord; CCoif[1] = 1/(2*critCoord); CConst[1] = flatChance;
        }
    }
    //printf("%f %f %f, %f %f, %f %f\n", CDiff[0], CDiff[1], CDiff[2], CCoif[0], CCoif[1], CConst[0], CConst[1]);
    // Initialize the random number generator
    srand(time(NULL));

    return loadn;
}

//         ****************PROGRAM CONTROL FUNCITONS***************

//Controls the physical simulation of a single particle
void Evolve ()
{   
    if (highp) {
        dist = 2;
        initParticleHighPressure(Particle);        
    } 
    else{
        initParticleLow (Particle);            //** Initialize Particle State
    }
	moveParticle (Particle);            //** Move the particle one step to begin with.

    while(nmov <= 100 || (nmov <= 1e4 && Particle.vz < -.01)) //** Move the particle, but only for a maximum of 10000 steps
    {
        if (Particle.z > L-2 || ((Particle.vz > 0) && (Particle.z > (kill_height))))      //** If the particle reaches the top of the lattice or has upward movement and is above the surface,
        {                                                                   //** then the particle is considered to have flown off the surface.
            return;
        }

        // If the particle has entered an occupied site or just passed through a stable location (ie. a location with a nearest neighbor and with the
        // required number of total neighbors), then deposit the particle at the previous location. Note: In the case of side-wall depostion (non-direct hit),
        // the particle only has a sideDepProb chance of actually depositing.
        if (occupied(Particle.x,Particle.y,Particle.z) || ((RandomDouble(0,1) < sideDepProb) && valid(Particle.xl,Particle.yl,Particle.zl)))
        {
            if (RandomDouble(0,1) <= P[Particle.impact])      //** Decide if the particle will actually deposit verses be reemitted.
    	    {
                if (ETCH)
                {
                    //****Note: Etch functionality in this respect has not been tested or reviewed since major program alterations
                    //    and is not guarenteed to function correctly.
                    etch (Particle.x,Particle.y,Particle.z);   //** NOTE: this does not work with the side interaction
                    checklocal(Particle.x,Particle.y,Particle.z);
                }
                else
                {
                    int z;
                    if (SoS)   //** Drop the particle down to the top of the surface in solid-on-solid mode
                        z = h(Particle.xl, Particle.yl) + 1;
                    else
                        z = Particle.zl;

                    deposit (Particle.xl, Particle.yl, z);
                    if (!valid(Particle.xl, Particle.yl, z))            //** A direct hit deposition may lead to a particle that is not stable, and may need to be recovered.
                       recover(Particle.xl, Particle.yl, z, false);     //** However, as this is a new particle, there is no need to check for instabilities created by its recovery,
                                                                        //** and the recover function is told to ignore that step.
                }

                // This is where diffusion is initiated. Only particles within a small box (of size side*2+1) centered around the incident particle are considered.
                // For each incident particle, skips many diffusion attempts are made at random locations within the bounding box. As empty locations and invalid movement
                // directions are counted as attempts, the value of skips is usually in the area of a couple thousand. As such, diffusion is a very time consuming process
                // only recomended for smaller simulations. A value of skips = 0 will effectively disable diffusion. IMPORTANT NOTE: Though explicit diffusion may be disabled,
                // diffusion settings are still used in the recover function, so make sure you are aware of what those settings are.
//                if (SoS)
//                {
//                    // The value of skips is meant to signify diffusion strength in the ballistic deposition mode. For solid-on-solid, particles are only picked from from the top of
//                    // the surface and always move to unoccupied locations. Without a subsequent adjustment of the number of diffusion itterations, this would result in a drastically
//                    // increased diffusion rate. Calculations assuming an average of 8 out of 26 available movement locations per diffusable particle in ballistic aggregation have lead
//                    // to the below adjustment and should result in a roughly comparable diffusion rate.
//                    for (int m = 0; m < skips; ++m)
//                    {
//                        skip (RandomInt(Particle.xl-side, Particle.xl+side), RandomInt(Particle.yl-side, Particle.yl+side), 0);  //** The 3rd argument is irrelevant and so is just put to 0
//                    }
//                }
//                else
//                {
//                    for (int m = 0; m < skips; ++m)
//                    {
//                        skip (RandomInt(Particle.xl-side, Particle.xl+side), RandomInt(Particle.yl-side, Particle.yl+side), RandomInt(Particle.zl-side, Particle.zl+side));
//                    }
//                }


                for(int m=0;m<skips;m++){                    
                    boltzdiff(RandomInt(Particle.xl-side, Particle.xl+side), RandomInt(Particle.yl-side, Particle.yl+side), RandomInt(Particle.zl-side, Particle.zl+side));                    
                }

                return;
            }
    	    else
    	    {
                emitParticle (Particle);                            //** Reemit the particle if it is not going to be deposited.
                moveParticle (Particle);                            //** Move the particle one step away from the surface. To ensure the particle does not deposit again in the same location
                                                                    //** due to side-wall sticking, the particle must move a total of two steps before the next deposition check.
                if (occupied(Particle.x, Particle.y, Particle.z))   //** However, it is possible that the particle will immediately run into another particle. If this happens, the
                {                                                   //** previous movement step must be reset to ensure that the particle does not deposit into an already occupied space.
                    Particle.x = Particle.xl;
                    Particle.y = Particle.yl;
                    Particle.z = Particle.zl;
                    nmov = 1;
                }
            }
        }
        moveParticle (Particle);    //** Continue particle movement if it is not near the surface, has failed to interact with nearby particles, or has just been re-emitted.
    }
    return;
}

//This is the highest level function for a single simulation. It controls the progression of the simulation and the calculation/output of data.
void main_program ()
{
    char fileName[100];

    // Initialize the new simulation
    printf ("--Starting Simulation--\n");
    int startn = Initialize() + 1;

    // Create a discription file which will allow someone to later reference all details about the setup of the simulation.
    // All information is taken from the actual simulation parameters, so there is no possiblity of user error.
    printf ("creating discription file...\n");
    sprintf (fileName, "%s_Description.txt", baseName);
    ofstream fs;
    if (startn == 1)
    {
        fs.open(fileName, ios_base::trunc);
        fs << "Version 4.7\n\n";
    }
    else
    {
        fs.open(fileName, ios_base::ate);
        fs << "\n************************\nLOG: Continuing at n = %i with parameters,\n\n";
    }

    fs << "--------SIMULATION SIZE----------\n";
    fs << "Lattice width:  " << N << "  spaces\n";
    fs << "Lattice height:  " << L << "  spaces\n";
    fs << "Simulation time:  " << T << "  ticks,  " << (float)R/(N*N) << "  particles/area/tick\n";
    fs << "\n-------PHYSICAL BEHAVIOR---------\n";
    fs << "Simulation type:  "; if(SoS) fs << "Solid-on-Solid"; else fs << "Ballistic Aggregation"; fs << "\n";
    fs << "Deposition type:  "; if(!ETCH) fs << "Deposition"; else fs << "Etching"; fs << "\n";
    fs << "Side-wall deposition likelihood:  "; if(sideDepProb <= 0) fs << "0"; else if(sideDepProb >= 1) fs << "100"; else fs << sideDepProb*100; fs << "%\n";
    fs << "Neighbors required for stability:  " << validLim << "\n";
    fs << "Diffusion and Recovery:\n";
    fs << "  Diffusion type:  "; if(DiffType == 1) fs << "Boltzman Escape"; else if(DiffType == 2) fs << "Linear Escape";
          else if (DiffType == 3) fs << "Sighted Coord"; else if (DiffType == 4) fs << "Sighted Height"; else fs << "Error: none selected"; fs << "\n";
    fs << "  Diffusion attempts per incident:  " << skips; if (SoS) fs << " equivalent"; fs << "\n";
    fs << "  Selection range from incident:  " << side << "\n";
    fs << "  Coordination contribution from far particles:  " << farCoord << "\n";
    if (DiffType == 1)
    {
        fs << "  Energy per neighboring bond:  " << En << " eV\n";
        fs << "  Bond number modifier*:  " << Ea << "\n";
        fs << "  Substrate temperature:  " << Ts << " K\n";
    }
    else if (DiffType == 2)
    {
        fs << "  Coordination where movement is guarenteed:  " << coord100 << "\n";
        fs << "  Coordination where movement becomes impossible:  " << coord0 << "\n";
    }
    else if (DiffType == 3)
    {
        fs << "  Coordination gain needed for guarenteed movement:  " << critCoord << "\n";
        fs << "  Flat potential movement likelihood:  "; if(flatChance <= 0) fs << "0"; else if(flatChance >= 1) fs << "100"; else fs << flatChance*100; fs << "%\n";
    }
    else if (DiffType == 4)
    {
        fs << "  Height loss needed for guarenteed movement:  " << critCoord << "\n";
        fs << "  Flat height movement likelihood:  "; if(flatChance <= 0) fs << "0"; else if(flatChance >= 1) fs << "100"; else fs << flatChance*100; fs << "%\n";
    }
    fs << "Reemission:  ";
    if(P[0] < 1)
    {
        fs << "Active\n";
        fs << "  Sticking coefficients:  "; for(int i = 0; i < 10; ++i){fs << P[i]; if (P[i] < 1) fs << ", "; else break;} fs << "\n";
        fs << "  Reemission type:  "; if (emitDist == 1) fs << "Reflection\n"; else if (emitDist == 2) fs << "Cosine\n"; else fs << "Uniform\n";
    }
    else fs << "Inactive\n";
    fs << "\n----DISTRIBUTION AND ROTATION----\n";
    fs << "Flux Distribution:\n";
    fs << "  Type:  ";
    if(dist == 1){
        fs << "Uniform\n";
        fs << "  Angles of incidence:  From " << angleInner << " to " << angleOuter << " degrees off normal\n";}
    else
        fs << "Cosine\n";
    fs << "Substrate orientation:\n";
    fs << "  Oblique angle:  " << OblqAngle << "  degrees\n";
    fs << "  Tilt swing:  ";
    if(tiltswing) {
        fs << "Active\n";
        fs << "    Swing angle:  " << TiltAngle << "  degrees\n";
        fs << "    Speed:  " << PartPerTiltSwing << "  particles per period\n";
        fs << "    Period offset:  " << thetaphase << "  degrees\n";}
    else
        fs << "Inactive\n";
    fs << "  Rotation:  ";
    if(rotation) {
        fs << "Active\n";
        fs << "    Speed:  " << PartPerRot << "  particles per period\n";
        fs << "    Direction:  "; if(direct == 1) fs << "clockwise"; else fs << "counter-clockwise"; fs << "\n";}
    else
        fs << "Inactive\n";
    fs << "  Rotation swing:  ";
    if(rotswing) {
        fs << "Active\n";
        fs << "    Swing angle:  " << RotAngle << "  degrees\n";
        fs << "    Speed:  " << PartPerRotSwing << "  particles per period\n";
        fs << "    Period offset:  " << phiphase << "  degrees\n";}
    else
        fs << "Inactive\n";
    fs << "\n------------TEMPLATE-------------\n";
    fs << "Template:  ";
    if(UseTemplate == 1){
        fs << "Square\n";
        fs << "  Square heights:  " << Hsq << "  lattice spaces\n";
        fs << "  Square widths:  " << Lsq << "  lattice spaces\n";
        fs << "  Square separations:  " << Dsq << "  lattice spaces\n";}
    else if(UseTemplate == 2 || UseTemplate == 3){
        if(UseTemplate == 2) fs << "Flat connical\n"; else fs << "Arc connical\n";
        fs << "  Bottom radius:  " << rbottom << " lattice spaces\n";
        fs << "  Top radius:  " << rtop << " lattice spaces\n";
        fs << "  Cone height:  " << hFCone << " lattice spaces\n";
        fs << "  Peak to peak distance:  " << dPeakToPeak << " lattice spaces\n";}
    else if(UseTemplate == 4)
        fs << "One\n";
    else if(UseTemplate == 5){
        fs << "Colloid\n";
        fs << "  Radius:  " << Rcoll << " lattice spaces\n";}
    else if(UseTemplate == 6)
        fs << "Custom Template\n";
    else
        fs << "Flat\n";
    fs.close ();

    // Create the statistics file. This file will remain open throughout the simulation and will store information
    // such as the interface width and local slope in a space delineated format.
    printf ("creating statistics file...\n");
    sprintf (fileName, "%s_Statistics.txt", baseName);
    FILE * ofile;
    if (startn == 1)
    {
        ofile = fopen (fileName, "w");
        fprintf(ofile, "N height IW slope avgCoord rmsdCoord avgFinal rmsdFinal avgSkip avgRecov\n");
    }
    else
    {
        ofile = fopen (fileName, "a");
    }
    printf ("initialization complete ... \nstarting simulation...\n");

//** Start Deposition **
    int days = 0;       //** variables for measuring the simulation running time.
    int hours = 0;
    double mins = 0;

    char option;        //** variables for controling pause functions
    bool done = false;
    bool quit = false;

    int i, j, n;

    RandomInitialise(rand()%30000, rand()%30000);     //** Initialize the primary random number generator

    // Main loop of the simulation. Each iteration here represents one cycle of the simulation.
    for (n=startn; n <= T; ++n)
    {
        if ((int)n>10000){
            highp = true;
        }
        start = clock();
        numCoord = 0; avgCoord = 0; rmsCoord = 0;
        foundcount = 0; recovered = 0;
        skipcount = 0; avgFinal = 0; rmsFinal = 0;
        avgDown = 0; avgUp = 0;
        height = 0; width  = 0; slope  = 0;

        // Secondary loop of the simulation. All of the particles for the cycle are deposited here.
        for (i = 0;i < R;++i)
        {
            Evolve ();
            // Check to see if the surface is about to grow (or be etched) outside of the lattice and shift accordingly.
            if (ETCH == 0 && Max >= L - 8)
                ShiftDown ();
            else if (ETCH == 1 && Max < L - 8)
                ShiftUp ();
        }

        // From here on, the rest of the cycle is data collection and output.
        for (i = 0;i < N;++i)
        {
            for (j = 0;j < N;++j)
            {
                height += H[j*N+i]; //average height
            }
        }
        height = height/double(N*N);
        for (i = 0;i < N;++i)
        {
            for (j = 0;j < N;++j)
            {
                width += pow(H[j*N+i] - height, 2); //interface width
                slope += pow(double(h(i+1,j)-h(i-1,j))/2.0 , 2); //RMS surface slope
            }
        }
        width = sqrt(width/double(N*N));            //** Refers to the Interface Width (or RMSD) of the surface.
        slope = sqrt(slope/double(N*N));            //** Refers to the Local Slope of the surface
        height = height + double(h0);               //** Refers to the average heigh of the surface

        if (numCoord != 0)
        {
            avgCoord = avgCoord/numCoord;           //** Refers to the average coordination (or difference) seen by particles attempting to diffuse
            rmsCoord = sqrt(rmsCoord/numCoord);     //** Refers to the RMSD of the coordination seen by particles attempting to diffuse
            preCoord = avgCoord;                    //** For the sake of simplicity, the previous cycle's average value is used in calculating the RMSD
        }
        if(skipcount != 0)
        {
            avgFinal = avgFinal/skipcount;          //** Refers to the average coordination (or difference) seen by particles that successfully diffuse
            rmsFinal = sqrt(rmsFinal/skipcount);    //** Refers to the RMSD of the coordination seen by particles that successfully diffuse
            preFinal = avgFinal;                    //** For the sake of simplicity, the previous cycle's average value is used in calculating the RMSD
            avgUp = avgUp/skipcount;                //** Refers to the percentage of successful diffusions where the particle lost coordination (ie. went up in energy). Difftypes 3 and 4 only.
            avgDown = avgDown/skipcount;            //** Refers to the percentage of successful diffusions where the particle gained coordination (ie. went down in energy). Difftypes 3 and 4 only.
        }

        //Statistics file output
        fprintf (ofile, "%i %f %f %f %f %f %f %f %f %f\n", (int) n, height, width, slope, avgCoord, rmsCoord, avgFinal, rmsFinal, (float)skipcount/R, (float)recovered/R);
                                                    //** Skipcount/R is the average number of successful diffusions that take place per incident particle
                                                    //** Recovered/R is the average number of particles that were successfully recovered per incident particle
        //Command window output
        printf ("n = %i ... min = %i ... max = %i\n", (int) n, (int) Min, (int) Max);
        printf ("n = %i ... avg skips = %0.3f ... percent down = %0.2f ... percent up = %0.2f\n", (int) n, (float)skipcount/R, avgDown, avgUp);
        printf ("n = %i ... coordination = %0.3f +- %0.3f ... final coord = %0.3f +- %0.3f\n", (int) n, avgCoord, rmsCoord, avgFinal, rmsFinal);
        printf ("n = %i ... avg found = %0.3f ... avg recovered = %0.3f ... total deleted = %i\n\n", (int) n, (float)foundcount/R, (float)recovered/R,  (int)deleted);
                                                    //** Foundcount/R is the average number of particles that needed to be recovered per incident particles. If this does not match recovered/R then
                                                    //**    it means some particles were deleted this cycle.
                                                    //** Deleted is the total number of particles, over the entire simulation, that have been deleted.

        // This is where the lattice images are taken. This is not done every
        // cycle, but only when the cycle number n is in the SaveTime array.
        //if (height >= HeightSave[HeightNumber])
        if (SaveTime[SaveNumber] == n)
        {
            // If selected, save 3D view
            if (sl)
            {
                printf ("Saving 3D lattice...\n");
                sprintf    (fileName, "%s_n%i.l", baseName, int (n));
    		    saveLat    (fileName);
            }
            // If selected, find and output a slice along the y-axis
            if (sly)
            {
                printf   ("Saving y-axis slice...\n");
                findHcs();
                sprintf    (fileName, "%s_n%i.cs", baseName, int (n));
                saveCrossSec3D (fileName);
            }
            // If selected, find and output center slice along x-axis
            if (slx)
            {
                printf   ("Saving x-axis slice...\n");
                findHcs2();
                sprintf    (fileName, "%s_n%i.xcs", baseName, int (n));
                saveCrossSec3D2 (fileName);
            }
            // If selected, output the overhead view
            if (slz)
            {
                printf   ("Saving overhead view...\n");
                sprintf    (fileName, "%s_n%i.z", baseName, int (n));
                saveHeight (fileName);
            }
            // Save porosity data
            if (spd)
	    {
	        printf   ("Saving porosity data...\n");
                sprintf    (fileName, "%s_n%i.p", baseName, int (n));
                savePorosity (fileName);	      
	    }

            ++SaveNumber;
            ++HeightNumber;
            printf   ("Done... Resuming simulation...\n");
        }

        // Add the cycle time to the total simulation time.
        mins += ((double)(clock() - start))/CLOCKS_PER_SEC/60;
        while (mins >= 60)
        {
            hours += 1;
            mins -= 60;
        }
        while (hours >= 24)
        {
            days += 1;
            hours -= 24;
        }

        // Here is the simulation pause functionality. If a key was pressed during the past cycle, kbhit() will return true and the simulation will pause.
        if (false/*kbhit()*/)
        {
            done = false;       //** The simulation will only unpause if done is set to true
            do{
                cout << "\nSimulation paused. What would you like to do?\n(f)inish simulation, (s)ave and exit, (r)esume simulation\n";
                cin.sync();         //** Clear character buffer
                do{
                    option = cin.get();             //** cin.get will pause and wait for input if there is no character in the buffer
                }while(!(option=='f' || option=='s' || option=='r'));       //** Only accept a valid command
                if(option == 'f')
                {
                    cout << "\nFinish simulation. The simulation will not be resumable.\nAre you sure (y,n)?\n";
                    cin.sync();
                    do{
                        option = cin.get();
                    }while(!(option=='y' || option=='n'));
                    if (option == 'y')              //** Only proceed if the action is confirmed. Otherwise return to the main menu.
                    {
                        quit = true;                //** Setting quit to true will exit out of the main simulation loop prematurely. Always use this command instead of terminating the application
                                                    //** as doing so would not allow for the completion details to be put in the description file.
                        done = true;
                    }
                }
                else if (option == 's')
                {
                    cout << "\nSave and exit. State will be saved and simulation may be resumed later.\nAre you sure (y,n)?\n";
                    cin.sync();
                    do{
                        option = cin.get();
                    }while(!(option=='y' || option=='n'));
                    if (option == 'y')
                    {
                        cout << "\nSaving simulation...\n";
                        saveSim(n);
                        quit = true;
                        done = true;
                    }
                }
                else
                {
                    cout << "\nResuming simulation...\n";
                    done = true;
                }
            }while(!done);
        }
        if (quit)
            break;
    }

    fclose (ofile);

    // Add some final statistics to the bottom of the description file.
    sprintf (fileName, "%s_Description.txt", baseName);
    ofile = fopen (fileName, "a");
    fprintf (ofile, "----------------------------\nLOG: As of n = %i,\n", (int) n);
    fprintf (ofile, "  Simulation time since last save: %i days, %i hours, %i minutes\n", (int) days, (int) hours, (int)mins);
    fprintf (ofile, "  Total deleted particles: %i\n", (int) deleted);
    fprintf (ofile, "  Average density of film: %0.3f\n", (double(n*R-deleted))/((height-7)*N*N));    //** This is just a rough estimate of the film density.
    if (WentOutside)
        fprintf (ofile, "\nALERT: Surface data went outside the lattice structure for some\n       portion of the simulation since the last save. Accuracy of surface\n       behavior in those areas cannot be guaranteed.\n");
    fclose (ofile);

    // Memory clean-up
    for (i = 0;i < N*N;++i)
        delete [] Surface [i];
    delete [] Surface;
    delete [] H;
    for (i = 0; i < N; ++i)
    {
        delete [] Hcs[i];
        delete [] Hcs2[i];
    }
    delete [] Hcs;
    delete [] Hcs2;

    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//*****************************************************************************************************************
//  *** ** End Function Definitions ** ***

// Main function of the program. Under normal circumstances, this is the only part of the code you should ever change.
// Use this space to set parameter values and que up simulations. You can change any parameters, including lattice
// size, between runs, but make sure that each run has a unique baseName to avoid file overwrites.
// Note: Loops and If statements are particularly helpful here if you are going to be running several similar simulations.
// Other Note: If using the Dev-C++ environment, go to Tools -> Compiler Options -> Settings -> Optimization, and turn both
// optimization options on and to the highest level. This will observably increase the speed of the simulation.
int main ()
{
    // Size
    N = 512;                       //** Due to a stmCalc detail (error?), N must be a power of 2 for the PSD functionality to work.
    L = 512;                       //** By definition, L must be a multiple of 8
    int times512 = N/512;
    R = 62500*times512*times512; //65536; //** Particles/cycle. This and the line above let me just change N and have the simulation time-scale remain constant.
    T = 10000;                        //** Total length of simulation
    // Images
	sl = false;                     //** Enable/disable individual lattice images
	sly = false;
	slx = true;
	slz = true;
	spd = false;
    // Distribution/rotation/template
    additional_height = 0;
    dist = 4;                       //** 1 = uniform, 2 = cosine, 3 = reversed cosine 4 = cospower
    cospower = (0.0);
	SoS = false;
    angleInner = 0;
    angleOuter = 0;
	OblqAngle = 0;
    UseTemplate = 0;                //** 0=flat, 1-square, 2-flat connical, 3-arc connical, 4-one, 5-colloid, 6-custom (template.z)
    tiltswing = false;
	rotation = false;
	PartPerRot = (int)round(R*.01);
	rotswing = false;
	// Partcle properties
    validLim = 8;
    farCoord = 0;
    sideDepProb = 1;
    P[0] = 1; P[1] = 1;
    emitDist = 1;                   //** 1 = reflection, 2 = cosine, 3 = uniform
    // Diffusion
    skips = 100;
    DiffType = 1;                   //1 -> Boltzman Escape, 2 -> Linear Escape, 3 -> Sighted Coord, 4 -> Sighted Height
    En = 0.05;
    Ea = 0.08;
    Ts = 298;
    side = 5;
    // Isoflux Parameters
    highp = true;
    spawn_height=2;
    spawn_std = 1; //this parameter will cause a deviation + or - from spawn_height
    kill_height=200;
    killHDiff = 60;  

    // Que
    if (SoS)
    {
        if (dist == 1)
        {
            if (rotation)
                sprintf(baseName, "N %i L %i T %i D %.3f CycPerRot %.2f NoDiff SoS RotUnif", (int) N, (int) L, (int) T, (float)R/(N*N), (float)PartPerRot/R);
            else
            {
                if (angleInner == 0 && angleOuter == 0)
                {
                    if (skips == 0)
                        sprintf(baseName, "N %i L %i T %i D %.3f NoDiff SoS NormUnif", (int) N, (int) L, (int) T, (float)R/(N*N));
                    else
                        sprintf(baseName, "N %i L %i T %i D %.3f Skips %i SoS NormUnif", (int) N,  (int)L,  (int)T, (float)R/(N*N), skips);
                }
                else
                    sprintf(baseName, "N %i L %i T %i D %.3f Ang %.0f-%.0f NoDiff SoS AngleUnif", (int) N,  (int)L, (int) T, (float)R/(N*N), angleInner, angleOuter);
            }
        }
        else if (dist == 2)
            sprintf(baseName, "N %i L %i T %i D %.3f Stick %.2f NoDiff SoS Cos", (int) N, (int) L, (int) T, (float)R/(N*N), P[0]);
        else
            sprintf(baseName, "N %i L %i T %i D %.3f NoDiff SoS RevCos", (int) N, (int) L, (int) T, (float)R/(N*N));
    }
    else
    {
        if (dist == 1)
        {
            if (rotation)
                sprintf(baseName, "N %i L %i T %i D %.3f CycPerRot %.2f NoDiff BA RotUnif", (int) N, (int) L, (int) T, (float)R/(N*N), (float)PartPerRot/R);
            else
            {
                if (angleInner == 0 && angleOuter == 0)
                {
                    if (skips == 0)
                        sprintf(baseName, "N %i L %i T %i D %.3f NoDiff BA NormUnif", (int) N, (int) L, (int) T, (float)R/(N*N));
                    else
                        sprintf(baseName, "N %i L %i T %i D %.3f Skips %i BA NormUnif", (int) N, (int) L, (int) T, (float)R/(N*N), skips);
                }
                else
                    sprintf(baseName, "N %i L %i T %i D %.3f Ang %.0f-%.0f NoDiff BA AngleUnif", (int) N, (int) L, (int) T, (float)R/(N*N), angleInner, angleOuter);
            }
        }
        else if (dist == 2)
            sprintf(baseName, "N %i L %i T %i D %.3f Stick %.2f NoDiff BA Cos", (int) N, (int) L, (int) T, (float)R/(N*N), P[0]);
        else
            sprintf(baseName, "N %i L %i T %i D %.3f NoDiff BA RevCos", (int)N, (int) L, (int) T, (float)R/(N*N));
    }

    main_program();


	return 0;
}
