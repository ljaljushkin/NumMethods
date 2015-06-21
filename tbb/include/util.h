#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "mmio.h"
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "memory.h"

#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

using namespace std;
using namespace tbb;

const double ZERO_IN_CRS = 0.000001;

struct mtxMatrix {
    int N;
    int NZ;

    double *Value;
    int *Col;
    int *Row;

    int *RowIndex;
};

class Multiplicator
{
	mtxMatrix A, B;
	vector<int>* columns;
	vector<double>* values;
	int *row_index;
	public:

	Multiplicator(mtxMatrix& _A, mtxMatrix& _B, vector<int>* &_columns, 
		vector<double>* &_values, int *_row_index) :
		A(_A), B(_B), columns(_columns), values(_values), row_index(_row_index)
	{}

	void operator()(const blocked_range<int>& r) const {
		int begin = r.begin();
		int end = r.end();
		int N = A.N;

		int i, j, k;
		int *temp = new int[N];

		for (i = begin; i < end; i++) {
			// i-я строка матрицы A
			// Обнуляем массив указателей на элементы
			memset(temp, -1, N * sizeof(int));
			// Идем по ненулевым элементам строки и заполняем массив указателей
			int ind1 = A.RowIndex[i], ind2 = A.RowIndex[i + 1];
			for (j = ind1; j < ind2; j++) {
				int col = A.Col[j];
				temp[col] = j; // Значит, что a[i, НОМЕР] лежит 
				// в ячейке массива Value с номером temp[НОМЕР]
			}
			// Построен индекс строки i матрицы A
			// Теперь необходимо умножить ее на каждую из строк матрицы BT
			for (j = 0; j < N; j++) {
				// j-я строка матрицы B
				double sum = 0;
				int ind3 = B.RowIndex[j], ind4 = B.RowIndex[j + 1];
				// Все ненулевые элементы строки j матрицы B
				for (k = ind3; k < ind4; k++) {
					int bcol = B.Col[k];
					int aind = temp[bcol];
					if (aind != -1)
					sum += A.Value[aind] * B.Value[k];
				}
				if (fabs(sum) > ZERO_IN_CRS) {
					columns[i].push_back(j);
					values[i].push_back(sum);
					row_index[i]++;
				}
			}
	}
	delete [] temp;
	}
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
