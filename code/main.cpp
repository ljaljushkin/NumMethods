#include <stdio.h>
#include "sparse.h"
#include "util.h"

const double EPSILON = 0.000001;

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    printf("Invalid input parameters\n");
    return 1;
  }
  
  int N = atoi(argv[1]);
  int NZ = atoi(argv[2]);
  char *mtxFileName = NULL;
  char *vecFileName = NULL;
  if (argc > 3 && argc < 6)
  {
      mtxFileName = argv[3];
      vecFileName = argv[4];
  }

  if ((NZ > N) || (N <= 0) || (NZ <= 0))
  {
    printf("Incorrect arguments of main\n");
    return 1;
  }

  crsMatrix A;
  double *x, *b, *bm;
  double timeM=0, timeM1=0, diff=0;

  if (ReadMatrix(A, mtxFileName) != 0)
  {
    GenerateRegularCRS(1, N, NZ, A);
    WriteMatrix(A, "mtx.txt");
  }
  if (ReadVector(&x, N, vecFileName) != 0)
  {
    GenerateVector(2, N, &x);
    WriteVector(x, N, "vec.txt");
  }

  InitializeVector(N, &b);
  Multiplicate(A, x, b, timeM);

  InitializeVector(N, &bm);
  SparseMKLMult(A, x, bm, timeM1);

  CompareVectors(b, bm, N, diff);

  if (diff < EPSILON)
    printf("OK\n");
  else
    printf("not OK\n");
  printf("%d %d\n", A.N, A.NZ);
  printf("%.3f %.3f\n", timeM, timeM1);

  FreeMatrix(A);
  FreeVector(&b);
  FreeVector(&bm);
  FreeVector(&x);

  return 0;
}