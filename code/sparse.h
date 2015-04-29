#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <time.h>
#include <stdlib.h>

struct crsMatrix
{
  int N;  // ������ ������� (N x N)
  int NZ; // ���-�� ��������� ���������

  // ������ �������� (������ NZ)
  double* Value;
  // ������ ������� �������� (������ NZ)
  int* Col;
  // ������ �������� ����� (������ N + 1)
  int* RowIndex; 
};

// ������� ��������� ����������� ������� � �� ������ x
// ��������� ������������ � ������ b
int Multiplicate(crsMatrix A, double *x, double *b, double &time);

#endif // __SPARSE_H__