/*-------------------Computional Fluid Dynamics First Project-----------------------*/
/*---------Author: Mohammad Jamal Razeghi---------Student Number: 9665611105--------*/
/*---------Solving Advection Equation Using:------------------------------------------

--------------  (1)First Order Upwind Method  ----------------------------------------
--------------  (2)LAX Method  -------------------------------------------------------
--------------  (3)Leap Frog Method  -------------------------------------------------
--------------  (4)Lax_Wendrof Method  -----------------------------------------------
--------------  (5)Crank Nichlson Method  --------------------------------------------
--------------  (6)Second Order Upwind Method  ---------------------------------------
**************  (7)Rusanov Method  ***************************************************/


//  Domain of x:  0 <= x <= 1


/*--------------------------------------------------------------------------------------
|                                        Initial condition                             |
--------------------------------------------------------------------------------------*/
//                          u(x,0)=cos(4*pi*x) 


/*---------------------------------------------------------------------------------------
|                                        VARIABLES                                       |
-----------------------------------------------------------------------------------------|
T     =  Time                                                                            |
x     =  Distance                                                                        |
u     =  Velocity                                                                        |
C     =  Spreading Velocity of Wave                                                      |                        
a     =  Starting Point of the Grid                                                      |
b     =  Ending Point of the Grid                                                        |
co    =  Courant Number                                                                  |
omega =  Traslative Viscousity Coefficient                                               |
dx    =  Spacing of the Grid                                                             |
dt    =  Time Step                                                                       |
IM    =  Number of Grid Points                                                           |
n     =  Number of Time Steps                                                            |
----------------------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>


using namespace std;

int main ()
{
	// DEFINING VARIABLES

	const double  C   = 0.1;
	const int     T   = 5;
	const double  pi  = 3.1415926535897;
	double        a   = 0;
	double        b   = 1;
	double        dx  = 0.01;
	double        co  = 0.9;
	double        dt  = co * dx / C;
	double      omega = 4 * pow(co,2) - pow(co,4);           //to reduce dissipation Err
	//double      omega = ((4 * pow(co,2) + 1)*(4-pow(co,2)))/5;  //to reduce dispersion Err
	int            n  = floor(T / dt) + 1;
	int           IM  = floor((b-a) / dx) + 1;
	double*        x  = new double [IM];
	x [0]             = a;
	for (int j=1 ; j <= IM-1 ; j++)
		x [j]         = x [j-1] + dx;
	double*        u1 = new double [IM-1];
	double*        u2 = new double [IM-1];
	double**       u  = new double* [n];

	for (int i=1 ; i <= n ; i++)
		u [i-1] = new double [IM-1];

	/*---------------------------------------------*/

    //INITIAL VALUE

    for (int j=0 ; j <= IM-1 ; j++)
    {
        u [0][j] = cos(4*pi*x[j]);
	}

	/*----------------------------------------*/

	//SOLVING ADVECTION EQUATION USING RUSANOV METHOD


	for (int i=0 ; i <= n-2 ; i++)
	{
		// STEP ONE
	    for (int j=0 ; j <= IM-1 ; j++)
		{
			if (j == IM-1)
			{
				u1 [j] = 0.5 * ( u [i][1] + u [i][j] ) - (co / 3) * (u [i][1] - u [i][j]);
			}
			else
			{
				u1 [j] = 0.5 * ( u [i][j+1] + u [i][j] ) - (co / 3) * (u [i][j+1] - u [i][j]);
			}
		}
		// STEP TWO
		for (int j=0 ; j <= IM-1 ; j++)
		{
			if (j == 0)
			{
				u2 [j] = u [i][j] - (2 * co / 3) * (u1 [j] - u1 [IM-2]);
			}
			else
			{
				u2 [j] = u [i][j] - (2 * co / 3) * (u1 [j] - u1 [j-1]);
			}
		}
		//STEP THREE
		for (int j=0 ; j <= IM-1 ; j++)
		{
			if (j == IM-1)
			{
				u [i+1][j] = u [i][j] - (co / 24) * (-2 * u [i][2] + 7 * u [i][1] - 7 * u [i][j-1] + 2 * u [i][j-2]) - (3 * co / 8) * (u2 [1] - u2 [j-1])
								- (omega / 24) * (u [i][2] - 4 * u [i][1] + 6 * u [i][j] - 4 * u [i][j-1] + u [i][j-2]);
			}
			else if (j == IM-2)
			{
				u [i+1][j] = u [i][j] - (co / 24) * (-2 * u [i][1] + 7 * u [i][j+1] - 7 * u [i][j-1] + 2 * u [i][j-2]) - (3 * co / 8) * (u2 [j+1] - u2 [j-1])
								- (omega / 24) * (u [i][1] - 4 * u [i][j+1] + 6 * u [i][j] - 4 * u [i][j-1] + u [i][j-2]);
			}
			else if (j == 0)
			{
				u [i+1][j] = u [i][j] - (co / 24) * (-2 * u [i][j+2] + 7 * u [i][j+1] - 7 * u [i][IM-2] + 2 * u [i][IM-3]) - (3 * co / 8) * (u2 [j+1] - u2 [IM-2])
								- (omega / 24) * (u [i][j+2] - 4 * u [i][j+1] + 6 * u [i][j] - 4 * u [i][IM-2] + u [i][IM-3]);
			}
			else if (j == 1)
			{
				u [i+1][j] = u [i][j] - (co / 24) * (-2 * u [i][j+2] + 7 * u [i][j+1] - 7 * u [i][j-1] + 2 * u [i][IM-2]) - (3 * co / 8) * (u2 [j+1] - u2 [j-1])
								- (omega / 24) * (u [i][j+2] - 4 * u [i][j+1] + 6 * u [i][j] - 4 * u [i][j-1] + u [i][IM-2]);
			}
			else
			{
				u [i+1][j] = u [i][j] - (co / 24) * (-2 * u [i][j+2] + 7 * u [i][j+1] - 7 * u [i][j-1] + 2 * u [i][j-2]) - (3 * co / 8) * (u2 [j+1] - u2 [j-1])
								- (omega / 24) * (u [i][j+2] - 4 * u [i][j+1] + 6 * u [i][j] - 4 * u [i][j-1] + u [i][j-2]);
			}
		}
	}

	/*------------------------------------------------------------------------------*/

	//CREATING TXT FILE OF OUTPUTS

	ofstream output;
		output.open("test.txt",ios::out);


		for (int i=0 ; i <= n-1 ; i++)
		{
		  for (int j=0 ; j <= IM-1 ; j++)
			{
				output << "\t" << u [i][j];
			}
		output << endl;
		}
	output.close();

	/*--------------------------------------------------------*/
	//Creating Tecplot file

	ofstream plot;
	plot.open("Rusanov2.txt",ios::out);
	//plot << "ZONE T=\"dx=" << dx << "\"" << endl;
	plot << "ZONE T=\"Co=" << co << "\"" << endl;
	for (int j=0 ; j <= IM-1 ; j++)
	{
		plot << x [j] << "\t" << u [n-1][j] <<endl;
	}
	cout << "omega = " << omega << endl;
	plot.close();
	delete [] u,x,u1,u2;
	return 0;
}
