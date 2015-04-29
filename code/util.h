#ifndef __UTIL_H__
#define __UTIL_H__

#include <memory.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "mkl_spblas.h"
#include "sparse.h"
#include <stdio.h>

// ������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// �������� ������ ��� ���� Value, Col � RowIndex
// ���������� ����� �������� mtx
void InitializeMatrix(int N, int NZ, crsMatrix &mtx);

// ������� ������ �� N ���������
// ���������� ����� �������� vec
void InitializeVector(int N, double **vec);

// ����������� ������ ��� ����� mtx
void FreeMatrix(crsMatrix &mtx);

// ����������� ������ ��� �������
void FreeVector(double **vec);

// ������� ����� imtx � omtx, ������� ������ ��� ���� Value, Col � RowIndex
void CopyMatrix(crsMatrix imtx, crsMatrix &omtx);

int CompareVectors(double* vec1, double* vec2, int n, double &diff);

int SparseMKLMult(crsMatrix A, double *x, double *b, double &time);

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// � ������ ������ cntInRow ��������� ���������
void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx);

// ���������� ������� ������
void GenerateVector(int seed, int N, double **vec);

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ����� ��������� ��������� � ������� ������ �� 1 �� cntInRow
// ����� ����� - ���������� ��������
void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx);

int WriteMatrix(crsMatrix mtx, char *fileName);

int ReadMatrix(crsMatrix &mtx, char *fileName);

int WriteVector(double *a, int N, char *fileName);

int ReadVector(double **a, int &N, char *fileName);

#endif // __UTIL_H__
