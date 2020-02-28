//New Boltzman Diffusion. Garauntees diffusion, if there is a possible location
int X2, Y2, Z2;
bool diff_flag;
void boltzdiff (int X1, int Y1, int Z1)
{
    
    if (!occupied(X1, Y1, Z1))
    { 
        return;
    }
    
    double E=RandomDouble(0,1); //particle "energy"
    double C=0;
    
    C=coordination(X1,Y1,Z1);
    
    if (!(E < exp(-(En*((double)C) + Ea))/(kB*Ts)))
    {
        return;        
    }
    
    int npaths=0;
    int paths[26];
    int path_RN=0;
    int path=0;
    
    unsigned char mask;
    int origH;

    //remove particle from surface

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
    
    
    for (int n = 0; n < 26; n++) 
    {
        X2 = X1 + S[n][0];
        Y2 = Y1 + S[n][1];
        Z2 = Z1 + S[n][2];
        if((!occupied(X2,Y2,Z2)) && (valid(X2,Y2,Z2)))
        {
            paths[npaths]=n;
            npaths++;
        }
    }
    
    if (npaths==0)
    {
        
        if (Z1 >= 0)             //**  return the particle to its position.
            Surface[Y1*N + X1][Z1/8] = Surface[Y1*N + X1][Z1/8] ^ mask;
        else
            H[Y1*N + X1] = origH;
        return;
	}
	
	else{
        
        path_RN=RandomInt(0,npaths-1);
        path = paths[path_RN];
        
        X2=X1+S[path][0]; 
        Y2=Y1+S[path][1]; 
        Z2=Z1+S[path][2];
        
        deposit(X2,Y2,Z2);
        checklocal(X1, Y1, Z1);     //** Since a particle left this location, check to ensure the                                    //** stability of neighbors.
        return;
	}
}


double mcoord(int x, int y, int z){
    double c=0.0;
  //matts model
    for (int i=6; i<=13;i++){
       if(h(x+S[i][0],y+S[i][1])>=z){
            c++;
       }
   }
   return c;
}

void boltzdiffSoS (int X1, int Y1, int Z1){
    
    while (X1 < 0)  X1 = X1 + N;
    while (Y1 < 0)  Y1 = Y1 + N;
    while (X1 >= N) X1 = X1 - N;
    while (Y1 >= N) Y1 = Y1 - N;
    
	if (!occupied(X1, Y1, Z1)){   //** Make sure that there is actually a particle to move
		return;
	}
    double E=RandomDouble(0,1); //particle "energy"
    double C=0;
    C=mcoord(X1,Y1,Z1);
    if (!(E <= exp(-(Ea+En*C)/(kB*Ts)))){
        return;
	}
	int npaths=0;
	int paths[8];
	int path_RN=0;
	int path_sum=0;
	int path=0;
	unsigned char mask;
	int origH;
	for (int m=0;m<8;m++){
		paths[m]=0;
	}
	//remove particle from surface


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
	///////////////////////////////
	for (int n=6;n<=13;n++){
		X2 = X1 + S[n][0];
		Y2 = Y1 + S[n][1];
		Z2 = h(X2,Y2)+1;
		if((Z2<=Z1) && valid(X2,Y2,Z2) && (!occupied(X2,Y2,Z2))){
		    //particles can diffuse at most 1 height unit above previous location, but can drop as far as possible
			//check path validity
            paths[npaths]=n;
			npaths++;
		}
	}
	if (npaths==0){
		if (Z1 >= 0)             //**  return the particle to its position.
			Surface[Y1*N + X1][Z1/8] = Surface[Y1*N + X1][Z1/8] ^ mask;
		else
			H[Y1*N + X1] = origH;
		return; // And leave the function
	}
	else{
		path_RN=RandomInt(0,npaths-1);
		path = paths[path_RN];
        
		X2=X1+S[path][0]; 
        Y2=Y1+S[path][1]; 
        Z2=h(X2,Y2)+1;
        
		etch(X1,Y1,Z1);
		deposit(X2,Y2,Z2);
		checklocal(X1, Y1, Z1);     //** Since a particle left this location, check to ensure the                                    //** stability of neighbors.
        return;
        //boltzdiff(X2,Y2,Z2);
	}
}
