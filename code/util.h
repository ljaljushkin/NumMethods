#ifndef __UTIL_H__
#define __UTIL_H__

#include <memory.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "sparse.h"
#include <stdio.h>
#include "mmio.h"

// ������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// �������� ������ ��� ���� Value, Col � Row
// ���������� ����� �������� mtx
void InitializeMatrix(int N, int NZ, mtxMatrix &mtx);

// ������� ������ �� N ���������
// ���������� ����� �������� vec
void InitializeVector(int N, double **vec);

// ����������� ������ ��� ����� mtx
void FreeMatrix(mtxMatrix &mtx);

// ����������� ������ ��� �������
void FreeVector(double **vec);

// ������� ����� imtx � omtx, ������� ������ ��� ���� Value, Col � Row
void CopyMatrix(mtxMatrix imtx, mtxMatrix &omtx);

int CompareVectors(double* vec1, double* vec2, int n, double &diff);

int SparseMKLMult(mtxMatrix A, double *x, double *b, double &time);

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// � ������ ������ cntInRow ��������� ���������
void GenerateRegularMTX(int seed, int N, int cntInRow, mtxMatrix& mtx);

// ���������� ������� ������
void GenerateVector(int seed, int N, double **vec);

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ����� ��������� ��������� � ������� ������ �� 1 �� cntInRow
// ����� ����� - ���������� ��������
void GenerateSpecialCRS(int seed, int N, int cntInRow, mtxMatrix& mtx);

int  WriteMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode);

int ReadMatrix(mtxMatrix &mtx, FILE *fileName);

int WriteVector(double *a, int N, char *fileName);

int ReadVector(double **a, int &N, char *fileName);

#endif // __UTIL_H__
