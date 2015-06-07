 #include <iostream>
#include "timer.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void GetParameter(char* filename, int& nx, int& ny);
void Outputfile(double d, char* filename);

int main(int argc, char* argv[])
{
    int		nx,ny;
    GetParameter(argv[1], nx, ny);

    double	x0 = 0,
            xEnd = 1,
            x_step = (xEnd - x0) / nx, //  ��� �� ����� x
            y0 = 0,
            yEnd = 1,
            y_step = (yEnd - y0) / ny, // ��� �� ����� y
            **M,	// nx x ny
            **B,	// nx x ny
            eps = 1e-11,
            omega = 2.0 - 10 * x_step;

    // ������������� ������
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
        }
    }

    double dmax=0, temp, dm=0;
    int i,j,k=0;
    Timer timer;
    double x_st_2=x_step*x_step, y_st_2 = y_step*y_step;

    timer.start();
    do {
        k++;
        printf("k=%i dm=%.16f eps=%.10f dm>eps=%i\n",k,dmax,eps,dmax>eps);
        dmax = 0; // ������������ ��������� �������� u
        for ( i=1; i<nx; i++ )
            for ( j=1; j<ny; j++ )
            {
                temp = M[i][j];
                M[i][j] =  ( 1 - omega)*M[i][j] + omega * 0.5 * (	y_st_2*(M[i-1][j] +
                                    M[i+1][j]) +
                                    x_st_2*(M[i][j-1] +
                                    M[i][j+1]) -
                                    x_st_2*y_st_2* B[i][j]
                                    ) / (x_st_2 + y_st_2);
                dm = fabs(temp-M[i][j]);
                if ( dmax < dm ) dmax = dm;
            }
    } while ( dmax > eps && k < 10000);
    timer.stop();

    float timeInSeconds = timer.getElapsed();

    double d_method=dmax;
    dmax=0;
    for ( i=0; i<nx+1; i++ )
        for ( j=0; j<ny+1; j++ )
        {
            temp=sin(x0+i*x_step)*sin(y0+j*y_step);//ideal
            dm=fabs(temp-M[i][j]);
            if ( dmax < dm ) dmax = dm;
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
