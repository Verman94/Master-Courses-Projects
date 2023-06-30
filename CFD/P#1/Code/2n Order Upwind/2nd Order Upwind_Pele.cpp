/*-------------------Computional Fluid Dynamics First Project-----------------------*/
/*---------Author: Mohammad Jamal Razeghi---------Student Number: 9665611105--------*/
/*---------Solving Advection Equation Using:------------------------------------------

--------------  (1)First Order Upwind Method  ----------------------------------------
--------------  (2)LAX Method  -------------------------------------------------------
--------------  (3)Leap Frog Method  -------------------------------------------------
--------------  (4)Lax_Wendrof Method  -----------------------------------------------
--------------  (5)Crank Nichlson Method  --------------------------------------------
**************  (6)Second Order Upwind Method  ***************************************
--------------  (7)Rusanov Method   -------------------------------------------------*/


//  Domain of x:  0 <= x <= 1


/*--------------------------------------------------------------------------------------
|                                        Initial condition                             |
--------------------------------------------------------------------------------------*/
//                          u(x,0)=1.0           x <= 0.1 
//                          u(x,0)=0.0    0.1 <= x <= 0.3
//                          u(x,0)=1.0           x >= 0.3 


/*--------------------------------------------------------------------------------------
|                                        VARIABLES                                     |
----------------------------------------------------------------------------------------
T   =  Time                                                                            |
x   =  Distance                                                                        |
u   =  Velocity                                                                        |
C   =  Spreading Velocity of Wave                                                      |                        
a   =  Starting Point of the Grid                                                      |
b   =  Ending Point of the Grid                                                        |
co  =  Courant Number                                                                  |
dx  =  Spacing of the Grid                                                             |
dt  =  Time Step                                                                       |
IM  =  Number of Grid Points                                                           |
IM1 =  Grid Points Until x = 0.1                                                       |
IM2 =  Grid Points Until x = 0.3                                                       |
n   =  Number of Time Steps                                                            |
--------------------------------------------------------------------------------------*/

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
	double        a   = 0;
	double        b   = 1;
	double        dx  = 0.01;
	double        co  = 1.9;
	double        dt  = co * dx / C;
	int            n  = floor(T / dt) + 1;
	int           IM  = floor((b-a) / dx) + 1;
	int           IM1 = floor((0.1-a) / dx) + 1;
	int           IM2 = floor((0.3-a) / dx) + 1;
	double*        x  = new double [IM];
	x [0]             = a;
	for (int j=1 ; j <= IM-1 ; j++)
		x [j]         = x [j-1] + dx;
	double**   u_pred = new double* [n];
	for (int i=1 ; i <= n ; i++)
		u_pred [i-1] = new double [IM-1];
	double**       u  = new double* [n];
	for (int i=1 ; i <= n ; i++)
		u [i-1] = new double [IM-1];

	/*---------------------------------------------*/

    //INITIAL VALUE

    for (int j=0 ; j < IM1-1 ; j++)
    {
        u [0][j] = 0;
    }
    for (int j=IM1-1 ; j <= IM2-1 ; j++)
    {
        u [0][j] = 1.0;
    }
    for (int j=IM2 ; j < IM ; j++)
    {
        u [0][j] = 0;
    }

	/*----------------------------------------*/

	//SOLVING ADVECTION EQUATION USING SECOND ORDER UPWIND METHOD

	for (int i=0 ; i <= n-2 ; i++)
	{
	

		//Predict
	    for (int j=0 ; j <= IM-1 ; j++)
		{
			//Periodic Boundray Condition
			if (j == 0)
			{
				u_pred [i+1][j] = u [i][j] - co * (u [i][j] - u [i][IM-2]);
			}
			else
			{
				u_pred [i+1][j] = u [i][j] - co * (u [i][j] - u [i][j-1]);
			}
		}

		//Correct
		for (int j=0 ; j <= IM-1 ; j++)
		{
			//Periodic Boundray Condition
			if (j == 0)
			{
				u [i+1][j] = 0.5 * ( u [i][j] + u_pred [i+1][j] - co * (u_pred [i+1][j] - u_pred [i+1][IM-2]) - co * (u [i][j] - 2 * u [i][IM-2] + u [i][IM-3]));
			}
			else if (j == 1)
			{
				u [i+1][j] = 0.5 * ( u [i][j] + u_pred [i+1][j] - co * (u_pred [i+1][j] - u_pred [i+1][j-1]) - co * (u [i][j] - 2 * u [i][j-1] + u [i][IM-2]));
			}
			else
			{
				u [i+1][j] = 0.5 * ( u [i][j] + u_pred [i+1][j] - co * (u_pred [i+1][j] - u_pred [i+1][j-1]) - co * (u [i][j] - 2 * u [i][j-1] + u [i][j-2]));
			}
		}
	}

	/*----------------------------------------------------------------------*/

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
	plot.open("Upwind2.txt",ios::out);
	//plot << "ZONE T=\"dx=" << dx << "\"" << endl;
	plot << "ZONE T=\"Co=" << co << "\"" << endl;
	for (int j=0 ; j <= IM-1 ; j++)
	{
		plot << x [j] << "\t" << u [n-1][j] <<endl;
	}
	plot.close();
	delete [] u,x,u_pred;
	return 0;
}
