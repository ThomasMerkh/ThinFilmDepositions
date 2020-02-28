// vector2 and particle implementation ...

#include<cmath>
using namespace std;


struct vector2
{
   long int x;
   long int y;
   long int z;
};

const vector2 zeroVector = {0,0,0};

int operator == (const vector2 & L,const vector2 & R)
{
   return (L.x*R.y == L.y*R.x) && (L.x*R.z == L.z*R.x) && (L.y*R.z == L.z*R.y) && (L.x*L.y*L.z == R.x*R.y*R.z);
}

vector2 operator + (vector2 L,const vector2 & R)
{
   L.x = L.x + R.x;
   L.y = L.y + R.y;
   L.z = L.z + R.z;

   return L;
}

vector2 operator - (vector2 L,const vector2 & R)
{
   L.x = L.x - R.x;
   L.y = L.y - R.y;
   L.z = L.z - R.z;

   return L;
}

long int operator * (const vector2 & L,const vector2 & R)
{
   return L.x*R.x + L.y*R.y + L.z*R.z;
}

vector2 operator * (const long int & L,vector2 R)
{
   R.x = R.x*L;
   R.y = R.y*L;
   R.z = R.z*L;

   return R;
}

vector2 operator * (vector2 L,const long int & R)
{
   L.x = L.x*R;
   L.y = L.y*R;
   L.z = L.z*R;

   return L;
}

vector2 operator / (vector2 L,const long int & R)
{
   L.x = L.x/R;
   L.y = L.y/R;
   L.z = L.z/R;

   return L;
}

vector2 cross (vector2 & L,vector2 & R)
{
   vector2 Rv;

   Rv.x = L.y*R.z - L.z*R.y;
   Rv.y = L.z*R.x - L.x*R.z;
   Rv.z = L.x*R.y - L.y*R.x;

   return Rv;
}

inline vector2 reflected(vector2 normal,vector2 incident)
{
   vector2 R = normal;

   R = R*(2*(normal*incident));

   R = (normal*normal)*incident - R;

   return R;
}

float dot (vector2 L,vector2 R)
{
   return float(L*R)/sqrt(float((L*L)*(R*R)));
}

float norm (vector2 L)
{
   return sqrt (float(L*L));
}

vector2 Z (vector2 & normal,vector2 & ref)
{
   return normal;
}

vector2 X (vector2 & normal,vector2 & ref)
{
   return cross (normal,ref);
}

vector2 Y (vector2 & normal,vector2 & ref)
{
   vector2 temp = X(normal,ref);

   return cross (normal,temp);
}
/*
void clip (vector2 & vect)
{
   while (max(abs(vect.x),abs(vect.y),abs(vect.z)) > 512)
   {
      vect = vect/2;
   }

   return;
}
*/

struct particle
{
   int xn;           //
   int yn;           //Starting Particle positions
   int zn;           //

   long int x;            //
   long int y;            //  Current Particle Position ...
   long int z;       //

   int xl;           //
   int yl;           //  Previous Particle Position ...
   int zl;      //

   double vx;         //
   double vy;         //  Particle Velosity
   double vz;         //

   int impact;       //

   int dx;           //
   int dy;           //Facilitative Particle State ...
   long int dz;      //

   int mfp;           //
   bool etch;            //

};

void setPInternal (particle & P)
{
   P.xl = P.x;
   P.yl = P.y;
   P.zl = P.z;

   if (P.vx > 0)
      P.dx = 1;
   else if (P.vx < 0)
      P.dx = -1;
   else
      P.dx = 0;

   if (P.vy > 0)
      P.dy = 1;
   else if (P.vy < 0)
      P.dy = -1;
   else
      P.dy = 0;

   if (P.vz > 0)
      P.dz = 1;
   else if (P.vz < 0)
      P.dz = -1;
   else
      P.dz = 0;
}
