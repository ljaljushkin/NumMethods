#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <time.h>
#include <stdlib.h>

struct mtxMatrix
{
  int N;  // Размер матрицы (N x N)
  int NZ; // Кол-во ненулевых элементов

  // Массив значений (размер NZ)
  double* Value;
  // Массив номеров столбцов (размер NZ)
  int* Col;
  // Массив индексов строк (размер NZ)
  int* Row; 
};

// Функция умножения разреженной матрицы А на вектор x
// Результат записывается в вектор b
int Multiplicate(mtxMatrix A, double *x, double *b, double &time);

#endif // __SPARSE_H__