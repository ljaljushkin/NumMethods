#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <time.h>
#include <stdlib.h>

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
};

// ������� ��������� ����������� ������� � �� ������ x
// ��������� ������������ � ������ b
int Multiplicate(mtxMatrix A, double *x, double *b, double &time);

#endif // __SPARSE_H__