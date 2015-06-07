#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "timer.h"
#include <cilk/cilk.h> 
#include <cilk/cilk_api.h> 
#include <cilk/reducer_max.h>


double fM(double x, double y,double x_step, double y_step);
double fB(double x, double y);


cilk::reducer< cilk::op_max<double> > max_res;
cilk::reducer< cilk::op_max<double> > max_local;

int main(int argc, char* argv[])
{
	int		nx,ny;
	GetParameter(argv[1], nx, ny);

	double	x0 = 0, 
			xEnd = 1,	// чтобы убралось два периода 
			x_step = (xEnd - x0) / nx, //  шаг по сетке x
			y0 = 0,
			yEnd = 1,
			y_step = (yEnd - y0) / ny, // шаг по сетке y
			**M,	// nx ny
			**B,	// nx ny
			eps = 1e-11,
			omega = 2.0 - 10 * x_step;

	// инициализация матриц
	M = new double*[nx+1];
	for (int i=0; i < nx+1; i++) 
	{
		M[i] = new double[ny+1];
		for (int j=0; j < ny+1; j++)
		{
			M[i][j] = sin(x0+i*x_step)*sin(y0+j*y_step);
			if (i!=0 && i!=nx && j!=0 && j!=ny)
			{
				M[i][j] += -0.1 + double((double)rand() / RAND_MAX * 0.2);	
			}
		}
	}

	B = new double*[nx+1];
	for (int i=0; i < nx+1; i++) 
	{
		B[i] = new double[ny+1];
		for (int j=0; j < ny+1; j++)
		{
			B[i][j] =  -2 * sin(x0+i*x_step) * sin(y0+j*y_step);
			//B[i][j] = fB( x0+i*x_step, y0+j*y_step );
		}
	}
	//
	double dmax=0, *temp, *dm, *d;
	dm = new double[nx-1];
	d = new double[nx-1];
	temp = new double[nx-1];
	int i, *j,k=0;
	
	j = new int[nx-1];
	Timer timer;
	double x_st_2=x_step*x_step, y_st_2 = y_step*y_step;
	timer.start();

	int wawe=0;
	int kol_iter=0;
	do {
		max_res.set_value(0);
		kol_iter++;
		cilk_sync;
		//for ( wawe=1; wawe<nx-2; wawe++ ) // wawe 1:nx-3
		//printf("growing wawe\n");
		for ( wawe=1; wawe<nx+1-2; wawe++ ) // wawe 1:nx-3
		{
			dm[wawe] = 0;
			
			cilk_for ( i=1; i<wawe+1; i++ ) // i 1:wawe
			{
				j[i] = wawe + 1 - i; // j wawe:1
				//printf("i=%i,j=%i ", i,j);
				temp[i] = M[i][j[i]];
				M[i][j[i]] =  ( 1 - omega)*M[i][j[i]] + omega * 0.5 * (	y_st_2*(M[i-1][j[i]] + 
								M[i+1][j[i]]) +
								x_st_2*(M[i][j[i]-1] + 
								M[i][j[i]+1]) - 
								x_st_2*y_st_2* B[i][j[i]]
								) / (x_st_2 + y_st_2);
				d[i] = fabs(temp[i]-M[i][j[i]]);

				//printf("d[%i]=%f\n",i,d[i]);
				//if ( dm[wawe] < d ) dm[wawe] = d;
				
				//printf("I'm %i, d0: i=%i, j=%i, w=%i,dm=%.16f\n",omp_get_thread_num(),i,j,wawe,dm[wawe]);
			} // конец параллельной области
			
			cilk_sync;
			max_local.set_value(0);
			cilk_for ( i=1; i< wawe+1; i++ ) // max wawe=nx-2
			{	 
				max_local->calc_max(d[i]);	
			} // конец параллельной области 
			dm[wawe]=max_local.get_value();
			//printf("dm[%i]=%f\n\n",wawe,dm[wawe]);
			//system("pause");

		}
		cilk_sync;
		//system("pause");
		
		//printf("tuda: \n");
		//for (int l1=1; l1<nx-1; l1++) printf("dm[%i]=%f\n",l1,dm[l1]);
		//system("pause");
		// стоячие волны
		dm[wawe]=0;
		//printf("wawe=%i\n",wawe);
		//for ( k=0; k<ny-nx+1; k++ ) // k 0:ny-nx
		//printf("state wawe\n");
		for ( k=0; k<ny+1-nx-1+1; k++ ) // k 0:ny-nx
		{
			cilk_for( i=1; i<wawe+1; i++ ) // i 1:wawe
			{
				j[i] = wawe + k + 1 - i; // j wawe+k:k+1
				//printf("i=%i,j=%i ", i,j);
				temp[i] = M[i][j[i]];
				M[i][j[i]] =  ( 1 - omega)*M[i][j[i]] + omega * 0.5 * (	y_st_2*(M[i-1][j[i]] + 
								M[i+1][j[i]]) +
								x_st_2*(M[i][j[i]-1] + 
								M[i][j[i]+1]) - 
								x_st_2*y_st_2* B[i][j[i]]
								) / (x_st_2 + y_st_2);
				d[i] = fabs(temp[i]-M[i][j[i]]);
				//if ( dm[wawe] < d ) dm[wawe] = d;
				//printf("st I'm %i, d0: i=%i, j=%i, w=%i, k=%i, dm=%.10f\n",omp_get_thread_num(),i,j,wawe,k,dm[wawe]);
			} // конец параллельной области
			max_local.set_value(0);
			cilk_for ( i=1; i< wawe+1; i++ ) // max wawe=nx-2
			{	 
				max_local->calc_max(d[i]);	
			} // конец параллельной области 
			dm[wawe]=max_local.get_value();
			//printf("\n");
		}
		cilk_sync;
		//printf("zdec: \n");
		//printf("dm[%i]=%f\n",wawe,dm[wawe]);
		//system("pause");
		// затухание волны
		//for ( wawe=nx-3; wawe>0; wawe-- ) // wawe nx-3:1
		//printf("lowing wawe\n");
		for ( wawe=nx+1-3; wawe>0; wawe-- ) // wawe nx-3:1
		{
			cilk_for( i=nx+1-1-wawe; i<nx+1-1; i++ ) // i nx-1-wawe:nx-2
			{ 
					//j = ny+nx-wawe-3-i; //  j ny-2:ny-1-wawe
				j[i] = ny+1+nx+1-wawe-3-i; //  j ny-2:ny-1-wawe
				//printf("i=%i,j=%i ", i,j);
				temp[i] = M[i][j[i]];
				M[i][j[i]] =  ( 1 - omega)*M[i][j[i]] + omega * 0.5 * (	y_st_2*(M[i-1][j[i]] + 
								M[i+1][j[i]]) +
								x_st_2*(M[i][j[i]-1] + 
								M[i][j[i]+1]) - 
								x_st_2*y_st_2* B[i][j[i]]
								) / (x_st_2 + y_st_2);
				d[i] = fabs(temp[i]-M[i][j[i]]);
				//if ( dm[wawe] < d ) dm[wawe] = d;
					
				//printf("I'm %i, d0: i=%i, j=%i, w=%i\n",omp_get_thread_num(),i,j,wawe);
			} // конец параллельной области 
			max_local.set_value(0);
			cilk_for ( i=nx+1-1-wawe; i< nx+1-1; i++ ) // max wawe=nx-2
			{	 
				max_local->calc_max(d[i]);	
			} // конец параллельной области 
			dm[wawe]=max_local.get_value();
			//printf("\n");
		} 
		
		//printf("obratno: \n");
		//for (int l1=nx-2; l1>0; l1--) printf("dm[%i]=%f\n",l1,dm[l1]);


		//system("pause");
		
			//for ( wawe=1; wawe<nx-1; wawe++ ) // max wawe=nx-2
		cilk_sync;
		cilk_for( wawe=1; wawe<nx+1-1; wawe++ ) // max wawe=nx-2
		{	 
			max_res->calc_max(dm[wawe]);	
		} // конец параллельной области 
		cilk_sync;	
		printf("k=%i dmax=%.10f\n",kol_iter,max_res.get_value());
	} while (  max_res.get_value() > eps );
	timer.stop();
	float timeInSeconds = timer.getElapsed();
	
	double d_cur,d_method= max_res.get_value();
	dmax=0;
	double temp1;
	for ( i=0; i<nx+1; i++ )
		for (int j=0; j<ny+1; j++ )
		{
			temp1=sin(x0+i*x_step)*sin(y0+j*y_step);//ideal
			d_cur=fabs(temp1-M[i][j]);
			if ( dmax < d_cur ) dmax = d_cur;
		}
	
	char* outputfile = argv[2];
	Outputfile(dmax , outputfile);

	double duration = timer.getElapsed();
	
	char* timefile = argv[3];
	double duration = timer.getElapsed();
	Outputfile(duration, timefile);

	delete[]M;
	delete M;
	delete[]B;
	delete B;

	return 0;
}

void Outputfile(double d, char* filename)
{
	FILE* ptFile = fopen(filename, "w");
	fprintf(ptFile, "%.30f", d);
	fclose(ptFile);
}

void GetParameter(char* filename, int& nx, int& ny)
{
	nx = 0;
	ny = 0;
	FILE* ptFile = fopen(filename, "r");
	int m = fscanf(ptFile, "%i %i", &nx, &ny);
	fclose(ptFile);
}
