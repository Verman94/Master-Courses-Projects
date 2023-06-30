/*-------------------Computional Fluid Dynamics First Project-----------------------*/
/*---------Author: Mohammad Jamal Razeghi---------Student Number: 9665611105--------*/
/*---------Solving Advection Equation Using:------------------------------------------

--------------  (1)First Order Upwind Method  ----------------------------------------
--------------  (2)LAX Method  -------------------------------------------------------
--------------  (3)Leap Frog Method  -------------------------------------------------
**************  (4)Lax_Wendrof Method  ***********************************************
--------------  (5)Crank Nichlson Method  --------------------------------------------
--------------  (6)Second Order Upwind Method  ---------------------------------------
--------------  (7)Rusanov Method   -------------------------------------------------*/


//  Domain of x:  0 <= x <= 1


/*--------------------------------------------------------------------------------------
|                                        Initial condition                             |
--------------------------------------------------------------------------------------*/
//                          u(x,0)=cos(4*pi*x)


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
	const double  pi  = 3.1415926535897;
	double        a   = 0;
	double        b   = 1;
	double        dx  = 0.01;
	double        co  = 1;
	double        dt  = co * dx / C;
	int            n  = floor(T / dt) + 1;
	int           IM  = floor((b-a) / dx) + 1;
	double*        x  = new double [IM];
	x [0]             = a;
	for (int j=1 ; j <= IM-1 ; j++)
		x [j]         = x [j-1] + dx;
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

	//SOLVING ADVECTION EQUATION USING Lax_Wendrof METHOD

	for (int i=0 ; i <= n-2 ; i++)
	{
	    for (int j=0 ; j <= IM-1 ; j++)
		{
			//PERIODIC BOUNDRAY CONDITION

			if (j == 0)
			{
				u [i+1][j] = u [i][j] - (co / 2) * (u [i][j+1] - u [i][IM-2]) + (pow(co,2) / 2) * (u [i][j+1] - 2 * u [i][j] + u [i][IM-2]);
			}

			else if (j == IM-1)

			{
				u [i+1][j] = u [i][j] - (co / 2) * (u [i][1] - u [i][j-1]) + (pow(co,2) / 2) * (u [i][1] - 2 * u [i][j] + u [i][j-1]);
			}

			else
			{
				u [i+1][j] = u [i][j] - (co / 2) * (u [i][j+1] - u [i][j-1]) + (pow(co,2) / 2) * (u [i][j+1] - 2 * u [i][j] + u [i][j-1]);
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
	plot.open("Lax_Wendrof2.txt",ios::out);
	//plot << "ZONE T=\"dx=" << dx << "\"" << endl;
	plot << "ZONE T=\"Co=" << co << "\"" << endl;
	for (int j=0 ; j <= IM-1 ; j++)
	{
		plot << x [j] << "\t" << u [n-1][j] <<endl;
	}
	plot.close();
	delete [] u,x;
	return 0;
}
