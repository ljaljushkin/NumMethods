#include "util.h"

using namespace std;
using namespace tbb;

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
		int Col = i; // ������� � AT - ������ � �
		for (j = j1; j < j2; j++) 
		{
			double V = imtx.Value[j];  // ��������
			int RIndex = imtx.Col[j];  // ������ � AT
			int IIndex = omtx.RowIndex[RIndex + 1];
			omtx.Value[IIndex] = V;
			omtx.Col  [IIndex] = Col;
			omtx.RowIndex[RIndex + 1]++;
		}
	}
}

// ��������� 2 ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ���������� C = A * B, C - � ������� CRS (3 �������, ���������� � ����)
//   ������ ��� C � ������ ��������� �� ����������
// ���������� ������� ���������� ��������: 0 - ��, 1 - �� ��������� ������� (N)
int Multiplicate(mtxMatrix A, mtxMatrix B, mtxMatrix &C)
{
	if (A.N != B.N)
	return 1;

	int N = A.N;
	int i;

	task_scheduler_init init();

	vector<int>* columns = new vector<int>[N];
	vector<double> *values = new vector<double>[N];
	int* row_index = new int[N + 1];
	memset(row_index, 0, sizeof(int) * N);

	int grainsize = 10;

	parallel_for(blocked_range<int>(0, A.N, grainsize),
		Multiplicator(A, B, columns, values, row_index));

	int NZ = 0;
	for(i = 0; i < N; i++) {
		int tmp = row_index[i];
		row_index[i] = NZ;
		NZ += tmp;
	}
	row_index[N] = NZ;

	InitializeMatrix(N, NZ, C);

	int count = 0;
	for (i = 0; i < N; i++) {
		int size = columns[i].size();
		memcpy(&C.Col[count], &columns[i][0], size * sizeof(int));
		memcpy(&C.Value[count], &values[i][0], size * sizeof(double));
		count += size;
	}
	memcpy(C.RowIndex, &row_index[0], (N + 1) * sizeof(int));

	delete [] row_index;
	delete [] columns;
	delete [] values;

	return 0;
}
