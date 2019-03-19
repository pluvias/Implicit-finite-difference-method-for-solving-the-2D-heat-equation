#include <cstdlib>
#include <iostream>
#include "math.h"
#include <fstream>

using std:: cin; using std:: cout; using std:: endl; using std::ofstream;



int main()
{
	//numerical method parameters
	
	const int Nx = 500; // nodes count in x-space
	const int Ny = 500; // nodes count in y-space
	double hx, hy, tau;
	hx = 1.0/Nx; // x variable step
	hy = 0.5/Ny; // y step
	
	tau = (hx > hy) ? hy :hx; // time step = min(hx,hy)
	
	const int Nt = 1000;

	// grid
	double *x = new double[Nx];  double *y = new double[Ny]; 
	x[0] = 0.0;  y[0] = 0.0;  

	for (int i = 0; (i < Nx) ; i++)
	{
		x[i + 1] = x[i] + hx;
	}
	for (int j=0; (j < Ny);j++)
	{
		y[j + 1] = y[j] + hy;
	}
	

	double **Up = new double*[Nx]; // storing values ​​from the previous step U(n+1/2)
	double **Up1 = new double*[Nx]; // U(n)
	double **Heat_Coeff = new double*[Nx]; // thermal conductivity (lambda)

	double A; // coeff-s of the tridiagonal matrix
	double C;   
	double B; 
	double F; // F - right side (system of linear equations)

	double **alpha = new double *[Nx]; // numerical method parameters for "grid function" computation
	double **betta = new double *[Nx]; 

	for(int i=0; i<Nx; i++)
	{
		alpha[i] = new double[Ny];    
		betta[i] = new double[Ny];

		Up[i] = new double[Ny];
		Up1[i] = new double[Ny];
		Heat_Coeff[i] = new double[Ny];
	}
	//==================================

	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			Up1[i][j]=300; // init. cond 
			Heat_Coeff[i][j] = ((x[i] >= 0.25) && (x[i] <= 0.65) && (y[j] >= 0.1) && (y[j] <=0.25)) ? 10e-2 : 10e-4;
			alpha[i][j] = 0.0;

			//Up[i][j] = Up1[i][j]; // copy for the numerical method
		}
	}

	//Upper and lower boundary conditions
	for (int j = 0; j < Ny; j++)
	{
		Up1[0][j] = 600;
		Up1[Nx-1][j] = 1200;
		//Up[0][j] = Up1[0][j];
		//Up[Nx-1][j] = Up1[Nx-1][j];
		betta[0][j] = Up1[0][j];
		betta[Nx-1][j] = Up1[Nx-1][j];
	}

	//left and right boundary conditions
	for (int i = 0; i < Nx; i++)
	{
		Up1[i][0] = 600 * (1 + x[i]);
		Up1[i][Ny-1]= 600 * (1+pow(x[i],3));
		//Up[i][0] = Up1[i][0];
		//Up[i][Ny-1] = Up1[i][Ny-1];
		betta[i][0] = Up1[i][0];
		betta[i][Ny-1] = Up1[i][Ny-1];		
	}

	

	//ofstream check;
	//check.open("lambda.txt");
	//
	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		check <<  Heat_Coeff[i][j] << "      ";
	//	}
	//	check << ";" << endl;
	//}
	//check.close();
	

	for(int time=0; time<Nt; time++) // time loop
	{
		cout << "time: " << time << endl << endl; 
			  
		for(int i=0; i< Nx; i++)
		{
			for(int j=0; j<Ny; j++)
			{
				Up[i][j] = Up1[i][j];
			}
		}


		for(int j=1; j< Ny-1; j++) // y-loop
		{	  
			for(int i=1; i< Nx-1; i++)
			{
				A = - (Heat_Coeff[i-1][j] + Heat_Coeff[i][j]) / (4*hx*hx);
				B = - (Heat_Coeff[i+1][j] + Heat_Coeff[i][j]) / (4*hx*hx);		
				C = 1/tau - A - B;	
				F = Up[i][j]/tau + ((Heat_Coeff[i][j+1] + Heat_Coeff[i][j]) *(Up[i][j+1]-Up[i][j])-(Heat_Coeff[i][j-1]+Heat_Coeff[i][j])*(Up[i][j]-Up[i][j-1]))/(4*hy*hy);
					
				alpha[i][j] = -B / (C+A*alpha[i-1][j]);
				betta[i][j] = (F - A*betta[i-1][j])/(C + A*alpha[i-1][j]);
				//cout << "Received betta: " << betta[i][j] << endl << endl; //getch();
			}
	
			for(int i=Nx-1; i>1; i--)
			{
				Up1[i-1][j]=alpha[i-1][j]*Up1[i][j]+betta[i-1][j]; // find Up
				//Up[i][j] = Up1[i][j];
			}
		}

		
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				Up[i][j] = Up1[i][j]; // remember array
			}	
		}
				
		for(int i=1; i<Nx-1; i++) // x-loop
		{
			for(int j=1; j< Ny-1; j++)
			{
				A = - (Heat_Coeff[i][j-1] + Heat_Coeff[i][j]) / (4*hy*hy);
				B = - (Heat_Coeff[i][j+1] + Heat_Coeff[i][j]) / (4*hy*hy);
				C = 1/tau - A - B;
				
				F = Up[i][j]/tau + ((Heat_Coeff[i+1][j] + Heat_Coeff[i][j])*(Up[i+1][j]-Up[i][j])-(Heat_Coeff[i-1][j]+Heat_Coeff[i][j])*(Up[i][j]-Up[i-1][j]))/(4*hx*hx);
				
				alpha[i][j] = -B / (C+A*alpha[i][j-1]);
				betta[i][j] = (F - A*betta[i][j-1])/(C + A*alpha[i][j-1]);
			}

			for(int j=Ny-1; j>1; j--)
			{
				Up1[i][j-1]=alpha[i][j-1]*Up1[i][j]+betta[i][j-1]; // find Up1
				//Up[i][j] = Up1[i][j];
			}
		}

	}

	ofstream outt;
	//outt.precision(3);
	outt.open("heat.txt");
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			outt <<  Up1[i][j] << "  ";
		}
		outt << ";" << endl;
	}
	outt.close();
	

	//system ("Pause");
	return 0;
}
