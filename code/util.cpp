#include "util.h"

const double MAX_VAL = 100.0;

// Создает квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Выделяет память под поля Value, Col и Row
// Возвращает через параметр mtx
void InitializeMatrix(int N, int NZ, mtxMatrix &mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ];
    mtx.Col = new int[NZ];
    mtx.Row = new int[NZ];
}

// Создает вектор из N элементов
// Возвращает через параметр vec
void InitializeVector(int N, double **vec)
{
    *vec = new double[N];
}

// Освобождает память, выделенную под поля mtx
void FreeMatrix(mtxMatrix &mtx)
{
    delete[] mtx.Value;
    delete[] mtx.Col;
    delete[] mtx.Row;
}

// Освобождает память для вектора
void FreeVector(double **vec)
{
    delete [](*vec);
}

// Создает копию imtx в omtx, выделяя память под поля Value, Col и Row
void CopyMatrix(mtxMatrix imtx, mtxMatrix &omtx)
{
  // Инициализация результирующей матрицы
  int N  = imtx.N;
  int NZ = imtx.NZ;
  InitializeMatrix(N, NZ, omtx);
  // Копирование
  memcpy(omtx.Value, imtx.Value, NZ * sizeof(double));
  memcpy(omtx.Col, imtx.Col, NZ * sizeof(int));
  memcpy(omtx.Row, imtx.Row, NZ * sizeof(int));
}

double next()
{
  return ((double)rand() / (double)RAND_MAX);
}

// Генерирует квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// В каждой строке cntInRow ненулевых элементов
void GenerateRegularMTX(int seed, int N, int cntInRow, mtxMatrix& mtx)
{
  int i, j, k, f, tmp, notNull, c;

  srand(seed);

  notNull = cntInRow * N;
  InitializeMatrix(N, notNull, mtx);

  for(i = 0; i < N; i++)
  {
    // Формируем номера столбцов в строке i
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
    // Сортируем номера столцов в строке i
    for (j = 0; j < cntInRow - 1; j++)
      for (k = 0; k < cntInRow - 1; k++)
        if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1])
        {
          tmp = mtx.Col[i * cntInRow + k];
          mtx.Col[i * cntInRow + k] = mtx.Col[i * cntInRow + k + 1];
          mtx.Col[i * cntInRow + k + 1] = tmp;
        }
  }

  // Заполняем массив значений
  for (i = 0; i < cntInRow * N; i++)
    mtx.Value[i] = next() * MAX_VAL;

  // Заполняем массив индексов строк
  c = 0;
  for (i = 0; i < N; i++)
  {
      for (j = 0; j < cntInRow; j++)
      {
          mtx.Row[i*cntInRow + j] = i;
      }
  }
}

// Генерирует плотный вектор
void GenerateVector(int seed, int N, double **vec)
{
    srand(seed);
    InitializeVector(N, vec);
    for (int i = 0; i < N; i++)
    {
        (*vec)[i] = next() * MAX_VAL;
    }
}

// Генерирует квадратную матрицу в формате CRS (3 массива, индексация с нуля)
// Число ненулевых элементов в строках растет от 1 до cntInRow
// Закон роста - кубическая парабола
void GenerateSpecialCRS(int seed, int N, int cntInRow, mtxMatrix& mtx)
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
    mtx.Row[i] = columns[i].size();
    for (unsigned int j = 0; j < columns[i].size(); j++)
    {
      mtx.Col[count] = columns[i][j];
      mtx.Value[count] = next();
      count++;
    }
  }
  mtx.Row[N] = columns[N-1].size();

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

int WriteMatrix(mtxMatrix mtx, char *fileName)
{
    FILE *f = fopen(fileName, "w+");
    fprintf(f, "%i\n", mtx.NZ);
    fprintf(f, "%i\n", mtx.N);
    for (int i = 0; i < mtx.NZ; i++)
    {
        fprintf(f, "%lf;%i;%i\n", mtx.Value[i], mtx.Col[i], mtx.Row[i]);
    }
    fclose(f);
    return 0;
}

int ReadMatrix(mtxMatrix &mtx, char *fileName)
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
        fscanf(f, "%lf;%i;%i", &(mtx.Value[i]), &(mtx.Col[i]), &(mtx.Row[i]));
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
