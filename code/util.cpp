#include "util.h"

const double MAX_VAL = 100.0;

// ������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// �������� ������ ��� ���� Value, Col � RowIndex
// ���������� ����� �������� mtx
void InitializeMatrix(int N, int NZ, crsMatrix &mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ];
    mtx.Col = new int[NZ];
    mtx.RowIndex = new int[N + 1];
}

// ������� ������ �� N ���������
// ���������� ����� �������� vec
void InitializeVector(int N, double **vec)
{
    *vec = new double[N];
}

// ����������� ������, ���������� ��� ���� mtx
void FreeMatrix(crsMatrix &mtx)
{
    delete[] mtx.Value;
    delete[] mtx.Col;
    delete[] mtx.RowIndex;
}

// ����������� ������ ��� �������
void FreeVector(double **vec)
{
    delete [](*vec);
}

// ������� ����� imtx � omtx, ������� ������ ��� ���� Value, Col � RowIndex
void CopyMatrix(crsMatrix imtx, crsMatrix &omtx)
{
  // ������������� �������������� �������
  int N  = imtx.N;
  int NZ = imtx.NZ;
  InitializeMatrix(N, NZ, omtx);
  // �����������
  memcpy(omtx.Value   , imtx.Value   , NZ * sizeof(double));
  memcpy(omtx.Col     , imtx.Col     , NZ * sizeof(int));
  memcpy(omtx.RowIndex, imtx.RowIndex, (N + 1) * sizeof(int));
}


double next()
{
  return ((double)rand() / (double)RAND_MAX);
}

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// � ������ ������ cntInRow ��������� ���������
void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
  int i, j, k, f, tmp, notNull, c;

  srand(seed);

  notNull = cntInRow * N;
  InitializeMatrix(N, notNull, mtx);

  for(i = 0; i < N; i++)
  {
    // ��������� ������ �������� � ������ i
    for(j = 0; j < cntInRow; j++)
    {
      do
      {
        mtx.Col[i * cntInRow + j] = rand() % N;
        f = 0;
        for (k = 0; k < j; k++)
          if (mtx.Col[i * cntInRow + j] == mtx.Col[i * cntInRow + k])
            f = 1;
      } while (f == 1);
    }
    // ��������� ������ ������� � ������ i
    for (j = 0; j < cntInRow - 1; j++)
      for (k = 0; k < cntInRow - 1; k++)
        if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1])
        {
          tmp = mtx.Col[i * cntInRow + k];
          mtx.Col[i * cntInRow + k] = mtx.Col[i * cntInRow + k + 1];
          mtx.Col[i * cntInRow + k + 1] = tmp;
        }
  }

  // ��������� ������ ��������
  for (i = 0; i < cntInRow * N; i++)
    mtx.Value[i] = next() * MAX_VAL;

  // ��������� ������ �������� �����
  c = 0;
  for (i = 0; i <= N; i++)
  {
    mtx.RowIndex[i] = c;
    c += cntInRow;
  }
}

// ���������� ������� ������
void GenerateVector(int seed, int N, double **vec)
{
    srand(seed);
    InitializeVector(N, vec);
    for (int i = 0; i < N; i++)
    {
        (*vec)[i] = next() * MAX_VAL;
    }
}

// ���������� ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ����� ��������� ��������� � ������� ������ �� 1 �� cntInRow
// ����� ����� - ���������� ��������
void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
  srand(seed);
  double end = pow((double)cntInRow, 1.0 / 3.0);
  double step = end / N;

  std::vector<int>* columns = new std::vector<int>[N];
  int NZ = 0;

  for (int i = 0; i < N; i++)
  {
    int rowNZ = int(pow((double(i + 1) * step), 3) + 1);
    NZ += rowNZ;
    int num1 = (rowNZ - 1) / 2;
    int num2 = rowNZ - 1 - num1;

    if (rowNZ != 0)
    {
      if (i < num1)
      {
        num2 += num1 - i;
        num1 = i;
        for(int j = 0; j < i; j++)
          columns[i].push_back(j);
        columns[i].push_back(i);
        for(int j = 0; j < num2; j++)
          columns[i].push_back(i + 1 + j);
      }
      else
      {
        if (N - i - 1 < num2)
        {
          num1 += num2 - (N - 1 - i);
          num2 = N - i - 1;
        }
        for (int j = 0; j < num1; j++)
          columns[i].push_back(i - num1 + j);
        columns[i].push_back(i);
        for (int j = 0; j < num2; j++)
          columns[i].push_back(i + j + 1);
      }
    }
  }

  InitializeMatrix(N, NZ, mtx);

  int count = 0;
  int sum = 0;
  for (int i = 0; i < N; i++)
  {
    mtx.RowIndex[i] = sum;
    sum += columns[i].size();
    for (unsigned int j = 0; j < columns[i].size(); j++)
    {
      mtx.Col[count] = columns[i][j];
      mtx.Value[count] = next();
      count++;
    }
  }
  mtx.RowIndex[N] = sum;

  delete [] columns;
}


int CompareVectors(double* vec1, double* vec2, int n, double &diff)
{
    diff = DBL_MIN;
    for(int i=0; i<n; i++) 
    {
        double curr_diff = abs(vec1[i] - vec2[i]);
        if (diff > curr_diff)
            diff = curr_diff;
    }
    return 0;
}

int SparseMKLMult(crsMatrix A, double *x, double *b, double &time)
{
    // �������� ��������� ��� ������ ������� MKL
    // ��������������� ������� A � B � �������
    int i, j, n = A.N;
    for (i = 0; i < A.NZ; i++)
        A.Col[i]++;
    for (j = 0; j <= n; j++)
    {
        A.RowIndex[j]++;
    }

    // ������������ �������, ����������� C = op(A) * B
    char trans = 'N'; // ������� � ���, op(A) = A - �� ����� ��������������� A

    clock_t start = clock();
    mkl_dcsrgemv(&trans, &n, A.Value, A.RowIndex, A.Col, x, b);
    clock_t finish = clock();
    // �������� � ����������� ���� ������� A
    for (i = 0; i < A.NZ; i++)
        A.Col[i]--;
    for (j = 0; j <= A.N; j++)
    {
        A.RowIndex[j]--;
    }

    time = double(finish - start) / double(CLOCKS_PER_SEC);

    return 0;
}

int WriteMatrix(crsMatrix mtx, char *fileName)
{
    FILE *f = fopen(fileName, "w+");
    fprintf(f, "%i\n", mtx.NZ);
    fprintf(f, "%i\n", mtx.N);
    for (int i = 0; i < mtx.NZ; i++)
    {
        fprintf(f, "%lf;%i\n", mtx.Value[i], mtx.Col[i]);
    }
    for (int i = 0; i < mtx.N + 1; i++)
    {
        fprintf(f, "%i\n", mtx.RowIndex[i]);
    }
    fclose(f);
    return 0;
}

int ReadMatrix(crsMatrix &mtx, char *fileName)
{
    int N, NZ;
    if (fileName == NULL) return -1;
    FILE *f = fopen(fileName, "r");
    if (f == NULL) return -1;
    fscanf(f, "%d", &NZ);
    fscanf(f, "%d", &N);
    InitializeMatrix(N, NZ, mtx);
    for (int i = 0; i < NZ; i++)
    {
        fscanf(f, "%lf;%i", &(mtx.Value[i]), &(mtx.Col[i]));
    }
    for (int i = 0; i < N + 1; i++)
    {
        fscanf(f, "%d", &(mtx.RowIndex[i]));
    }
    fclose(f);
    return 0;
}

int WriteVector(double *a, int N, char *fileName)
{
    FILE *f = fopen(fileName, "w+");
    fprintf(f, "%i\n", N);
    for (int i = 0; i < N; i++)
    {
        fprintf(f, "%lf\n", a[i]);
    }
    fclose(f);
    return 0;
}

int ReadVector(double **a, int &N, char *fileName)
{
    if (fileName == NULL) return -1;
    FILE *f = fopen(fileName, "r");
    if (f == NULL) return -1;
    fscanf(f, "%i\n", &N);
    *a = new double[N];
    for (int i = 0; i < N; i++)
    {
        fscanf(f, "%lf\n", &((*a)[i]));
    }
    fclose(f);
    return 0;
}
