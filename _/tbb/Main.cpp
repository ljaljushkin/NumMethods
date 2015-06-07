#include "stdafx.h"

void GetParameter(char* filename, int& nx, int& ny);
void Outputfile(double d, char* filename);

int main(int argc, char* argv[])
{
	int		nx,ny;
	GetParameter(argv[1], nx, ny);

	task_scheduler_init init;
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
	
	//printf("done\n");
	//printf("allocating memory for B\n");
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

	double dmax=0, temp, *dm,d=0;
	//dm = new double[nx-2];
	dm = new double[nx+1-2];
	int i,j,k=0;
	Timer timer;
	double x_st_2=x_step*x_step, y_st_2 = y_step*y_step;
	timer.start();


	//omp_lock_t dmax_lock;
	//omp_init_lock(&dmax_lock);
	int wawe=0;
	int kol_iter=0;
	//CalcError calculating(dm,nx);
	CalcError calculating(dm,nx+1);
	do {
		dmax = 0; // максимальное изменение значений u
		// нарастание волны (nx – размер волны)
		kol_iter++;
		//for ( wawe=1; wawe<nx-2; wawe++ ) // wawe 1:nx-3
		for ( wawe=1; wawe<nx+1-2; wawe++ ) // wawe 1:nx-3
		{
			dm[wawe] = 0;
			/*parallel_for(	blocked_range<int>( 1, wawe+1, 40 ),	
							WawePropagation(M,B,wawe,dm,0,0,x_step,y_step,nx,ny) );*/
			parallel_for(	blocked_range<int>( 1, wawe+1, 40 ),	
							WawePropagation(M,B,wawe,dm,0,0,x_step,y_step,nx+1,ny+1,omega) );
		}
		// стоячие волны
		//system("pause");
		dm[wawe]=0;
		//for ( k=0; k<ny-nx+1; k++ ) // k 0:ny-nx
		for ( k=0; k<ny+1-nx-1+1; k++ ) // k 0:ny-nx
		{	
			/*parallel_for(	blocked_range<int>( 1, wawe+1, 40 ),	
							WawePropagation(M,B,wawe,dm,1,k,x_step,y_step,nx,ny) );*/
			parallel_for(	blocked_range<int>( 1, wawe+1, 40 ),	
							WawePropagation(M,B,wawe,dm,1,k,x_step,y_step,nx+1,ny+1,omega ) );
		}
		// затухание волны
		//system("pause");
		//for ( wawe=nx-3; wawe>0; wawe-- ) // wawe nx-3:1
		for ( wawe=nx+1-3; wawe>0; wawe-- ) // wawe nx-3:1
		{
			/*parallel_for(	blocked_range<int>( nx-1-wawe, nx-1, 40 ),	
							WawePropagation(M,B,wawe,dm,2,0,x_step,y_step,nx,ny) );*/
			parallel_for(	blocked_range<int>( nx+1-1-wawe, nx+1-1, 40 ),	
							WawePropagation(M,B,wawe,dm,2,0,x_step,y_step,nx+1,ny+1,omega ) );
		} 
		//system("pause");
		//CalcError calculating(dm,nx);
		CalcError calculating(dm,nx+1);
		//parallel_reduce(blocked_range<int>( 1, nx-1, 40 ), calculating);
		parallel_reduce(blocked_range<int>( 1, nx+1-1, 40 ), calculating);
		dmax = calculating.Result();
		
		printf("k=%i dmax=%.10f\n",kol_iter,dmax);
	} while ( dmax > eps );

	timer.stop();
	float timeInSeconds = timer.getElapsed();
	
	double dm1,d_method=dmax;
	dmax=0;
	for ( i=0; i<nx+1; i++ )
		for ( j=0; j<ny+1; j++ )
		{
			temp=sin(x0+i*x_step)*sin(y0+j*y_step);//ideal
			dm1=fabs(temp-M[i][j]);
			if ( dmax < dm1 ) dmax = dm1;
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
