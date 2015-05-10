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

// Создает квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Выделяет память под поля Value, Col и Row
// Возвращает через параметр mtx
void InitializeMatrix(int N, int NZ, mtxMatrix &mtx);

// Создает вектор из N элементов
// Возвращает через параметр vec
void InitializeVector(int N, double **vec);

// Освобождает память для полей mtx
void FreeMatrix(mtxMatrix &mtx);

// Освобождает память для вектора
void FreeVector(double **vec);

// Создает копию imtx в omtx, выделяя память под поля Value, Col и Row
void CopyMatrix(mtxMatrix imtx, mtxMatrix &omtx);

int CompareVectors(double* vec1, double* vec2, int n, double &diff);

int SparseMKLMult(mtxMatrix A, double *x, double *b, double &time);

// Генерирует квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// В каждой строке cntInRow ненулевых элементов
void GenerateRegularMTX(int seed, int N, int cntInRow, mtxMatrix& mtx);

// Генерирует плотный вектор
void GenerateVector(int seed, int N, double **vec);

// Генерирует квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Число ненулевых элементов в строках растет от 1 до cntInRow
// Закон роста - кубическая парабола
void GenerateSpecialCRS(int seed, int N, int cntInRow, mtxMatrix& mtx);

int  WriteMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode);

int ReadMatrix(mtxMatrix &mtx, FILE *fileName);

int WriteVector(double *a, int N, char *fileName);

int ReadVector(double **a, int &N, char *fileName);

#endif // __UTIL_H__
