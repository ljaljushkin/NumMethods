#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "mmio.h"

struct mtxMatrix {
    int N;
    int NZ;

    double *Value;
    int *Col;
    int *Row;

    int *RowIndex;
};

void InitializeMatrix(int N, int NZ, mtxMatrix &mtx);

void FreeMatrix(mtxMatrix &mtx);

void TriangleToFull(mtxMatrix *mtx, mtxMatrix *fullMatrix);

int WriteMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode);

int WriteFullMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode);

int ReadMatrix(mtxMatrix &mtx, FILE *fileName);

void Transpose(mtxMatrix imtx, mtxMatrix &omtx);

int Multiplicate(mtxMatrix A, mtxMatrix B, mtxMatrix &C);

#endif // __UTIL_H__
