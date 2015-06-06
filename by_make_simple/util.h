#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "mmio.h"

struct mtxMatrix
{
  int N;  // ������ ������� (N x N)
  int NZ; // ���-�� ��������� ���������

  // ������ �������� (������ NZ)
  double* Value;
  // ������ ������� �������� (������ NZ)
  int* Col;
  // ������ �������� ����� (������ NZ)
  int* Row; 

  int* RowIndex;
};

void InitializeMatrix(int N, int NZ, mtxMatrix &mtx);
void FreeMatrix(mtxMatrix &mtx);
int  WriteMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode);
int ReadMatrix(mtxMatrix &mtx, FILE *fileName);
#endif // __UTIL_H__
