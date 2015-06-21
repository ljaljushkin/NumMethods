#include <vector>
#include <stdlib.h>
#include <math.h>
#include "memory.h"

#include "util.h"

using namespace std;

const double ZERO_IN_CRS = 0.000001;

void InitializeMatrix(int N, int NZ, mtxMatrix &mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ];
    mtx.Col = new int[NZ];
    mtx.Row = new int[NZ];
    mtx.RowIndex = new int[N + 1];
}

void TriangleToFull(mtxMatrix* mtx, mtxMatrix* fullMatrix)
{
    int* insertionIndex = new int[fullMatrix->N];

    for(int i = 0; i < fullMatrix->N + 1; i++)
    {
         fullMatrix->RowIndex[i] = mtx->RowIndex[i];
    }

    for(int i = 0; i < fullMatrix->N; i++)
    {
        for(int j = mtx->RowIndex[i] + 1; j < mtx->RowIndex[i + 1]; j++)
        {
            for(int k = mtx->Col[j] + 1; k < fullMatrix->N + 1; k++)
            {
                fullMatrix->RowIndex[k]++;
            }
        }
    }

    for(int i = 0; i < fullMatrix->N; i++)
    {
         insertionIndex[i] = fullMatrix->RowIndex[i];
    }

    for(int i = 0; i < fullMatrix->N; i++)
    {
        for(int j = mtx->RowIndex[i]; j < mtx->RowIndex[i + 1]; j++)
        {
            fullMatrix->Col[insertionIndex[i]] = mtx->Col[j];
            fullMatrix->Value[insertionIndex[i]] = mtx->Value[j];
            insertionIndex[i]++;
            if (mtx->Col[j] != i)
            {
                fullMatrix->Col[insertionIndex[mtx->Col[j]]] = i;
                fullMatrix->Value[insertionIndex[mtx->Col[j]]] = mtx->Value[j];
                insertionIndex[mtx->Col[j]]++;
            }
        }
    }
}

void FreeMatrix(mtxMatrix &mtx)
{
    delete[] mtx.Value;
    delete[] mtx.Col;
    delete[] mtx.Row;
    delete[] mtx.RowIndex;
}

int WriteMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode)
{
    mm_write_banner(output_file, matcode);
    mm_write_mtx_crd_size(output_file, output_mtx.N, output_mtx.N, output_mtx.NZ);
    for (int i=0; i<output_mtx.NZ; i++)
        fprintf(output_file, "%d %d %.10lf\n", output_mtx.Row[i]+1, output_mtx.Col[i]+1, output_mtx.Value[i]);

    return 0;
}

int WriteFullMatrix(mtxMatrix &output_mtx, FILE *output_file, MM_typecode matcode)
{
    matcode[3] = 'G';
    mm_write_banner(output_file, matcode);
    mm_write_mtx_crd_size(output_file, output_mtx.N, output_mtx.N, output_mtx.NZ);
    for (int i = 0; i < output_mtx.N; i++)
        for(int j = output_mtx.RowIndex[i]; j < output_mtx.RowIndex[i + 1]; j++)
        {
            fprintf(output_file, "%d %d %.10lf\n", i + 1, output_mtx.Col[j] + 1, output_mtx.Value[j]);
        }

    return 0;
}

int ReadMatrix(mtxMatrix &mtx, FILE* inputfile) {
    int M, N, nz;
    if (mm_read_mtx_crd_size(inputfile, &M, &N, &nz) != 0)
        exit(1);

    InitializeMatrix(N, nz, mtx);

    for (int i = 0; i < nz; i++) {
        fscanf(inputfile, "%i %i %lf", &(mtx.Col[i]), &(mtx.Row[i]), &(mtx.Value[i]));
        mtx.Col[i]--;
        mtx.Row[i]--;
    }
    fclose(inputfile);
    return 0;
}

void Transpose(mtxMatrix imtx, mtxMatrix &omtx)
{
	int i, j;

	InitializeMatrix(imtx.N, imtx.NZ, omtx);

	memset(omtx.RowIndex, 0, (imtx.N + 1) * sizeof(int));
	for (i = 0; i < imtx.NZ; i++) 
		omtx.RowIndex[imtx.Col[i] + 1]++;
  
	int S = 0;
	for (i = 1; i <= imtx.N; i++) 
	{
		int tmp = omtx.RowIndex[i];
		omtx.RowIndex[i] = S;
		S = S + tmp;
	}

	for (i = 0; i < imtx.N; i++) 
	{
		int j1 = imtx.RowIndex[i];
		int j2 = imtx.RowIndex[i+1];
		int Col = i; // Столбец в AT - строка в А
		for (j = j1; j < j2; j++) 
		{
			double V = imtx.Value[j];  // Значение
			int RIndex = imtx.Col[j];  // Строка в AT
			int IIndex = omtx.RowIndex[RIndex + 1];
			omtx.Value[IIndex] = V;
			omtx.Col  [IIndex] = Col;
			omtx.RowIndex[RIndex + 1]++;
		}
	}
}

// Принимает 2 квадратных матрицы в формате CRS (3 массива, индексация с нуля)
// Возвращает C = A * B, C - в формате CRS (3 массива, индексация с нуля)
//   Память для C в начале считается не выделенной
// Возвращает признак успешности операции: 0 - ОК, 1 - не совпадают размеры (N)
int Multiplicate(mtxMatrix A, mtxMatrix B, mtxMatrix &C)
{
	if (A.N != B.N)
		return 1;

	int N = A.N;
	int i, j, k;

	vector<int> columns;
	vector<double> values;
	vector<int> row_index;

	int NZ = 0;
  
	int *temp = new int[N];
	row_index.push_back(0);

	for (i = 0; i < N; i++)
	{
		// i-я строка матрицы A
		// Обнуляем массив указателей на элементы
		memset(temp, -1, N * sizeof(int));
		// Идем по ненулевым элементам строки и заполняем массив указателей
		int ind1 = A.RowIndex[i], ind2 = A.RowIndex[i + 1];
		for (j = ind1; j < ind2; j++)
		{
			int col = A.Col[j];
			temp[col] = j; // Значит, что a[i, НОМЕР] лежит 
							// в ячейке массива Value с номером temp[НОМЕР]
		}
		// Построен индекс строки i матрицы A
		// Теперь необходимо умножить ее на каждую из строк матрицы BT
		for (j = 0; j < N; j++)
		{
			// j-я строка матрицы B
			double sum = 0;
			int ind3 = B.RowIndex[j], ind4 = B.RowIndex[j + 1];
			// Все ненулевые элементы строки j матрицы B
			for (k = ind3; k < ind4; k++)
			{
			int bcol = B.Col[k];
			int aind = temp[bcol];
			if (aind != -1)
				sum += A.Value[aind] * B.Value[k];
			}
			if (fabs(sum) > ZERO_IN_CRS)
			{
				columns.push_back(j);
				values.push_back(sum);
				NZ++;
			}
		}
		row_index.push_back(NZ);
	}

	InitializeMatrix(N, NZ, C);

	for (j = 0; j < NZ; j++)
	{
		C.Col[j] = columns[j];
		C.Value[j] = values[j];
	}
	for(i = 0; i <= N; i++)
		C.RowIndex[i] = row_index[i];

	delete [] temp;

	return 0;
}
