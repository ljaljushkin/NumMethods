#include <stdio.h>
#include "sparse.h"
#include "util.h"
#include "mmio.h"

const double EPSILON = 0.000001;


int ilu0(mtxMatrix &A, double * luval, int * uptr)
{
     int j1, j2; // граница текущей строки
     int jrow; // номер текущего столбца
     int k, j, jj; // счетчики циклов
     int *iw = NULL; // временный массив
     int jw;
     double t1;

     iw = new int[A.N];
     memset(iw, 0, A.N * sizeof(int));

     memcpy(luval, A.Value, A.RowIndex[A.N] * sizeof(double));

     for(k = 0; k < A.N; k++)
     {
        j1 = A.RowIndex[k];
        j2 = A.RowIndex[k + 1];
        for(j = j1; j < j2; j++)
        {
            iw[A.Col[j]] = j;
        }
		for(j = j1; (j < j2) && (A.Col[j] < k); j++)
		{
			jrow = A.Col[j];
			t1 = luval[j] / luval[uptr[jrow]];
			luval[j] = t1;

			for(jj = uptr[jrow]+1; jj < A.RowIndex[jrow + 1]; jj++)
			{
				jw = iw[A.Col[jj]];
				if(jw != 0)
				{
					luval[jw] = luval[jw] - t1 * luval[jj];
				}
			}
		}
		jrow = A.Col[j];
		uptr[k] = j;
		if((jrow != k) || (fabs(luval[j]) < EPSILON))
		{
			break;
		}
		for(j = j1; j < j2; j++)
		{
			iw[A.Col[j]] = 0;
		}
	 }

	 delete [] iw;
	 if(k < A.N)
		 return -(k+1);
	 return 0;
}

// swap entries in array v at positions i and j; used by quicksort
void swap(int * v, int i, int j)
{
  int t = v[i];
  v[i] = v[j];
  v[j] = t;
}

void swapd(double * v, int i, int j)
{
  double t = v[i];
  v[i] = v[j];
  v[j] = t;
}


// (quick) sort slice of array v; slice starts at s and is of length n
void quicksort(int * row, int * col, double * val, int s, int n)
{
  int x, p, i;
  // base case?
  if (n <= 1)
    return;
  // pick pivot and swap with first element
  x = row[s + n/2];
  swap(row, s, s + n/2);
  swap(col, s, s + n/2);
  swapd(val, s, s + n/2);
  // partition slice starting at s+1
  p = s;
  for (i = s+1; i < s+n; i++)
    if (row[i] < x) {
      p++;
      swap(row, i, p);
	  swap(col, i, p);
	  swapd(val, i, p);
    }
  // swap pivot into place
  swap(row, s, p);
  swap(col, s, p);
  swapd(val, s, p);
  // recurse into partition
  quicksort(row, col, val, s, p-s);
  quicksort(row, col, val, p+1, s+n-p-1);
}

void colQuickSort(int * col, double * val, int s, int n) {
  int x, p, i;
  // base case?
  if (n <= 1)
    return;
  // pick pivot and swap with first element
  x = col[s + n/2];
  swap(col, s, s + n/2);
  swapd(val, s, s + n/2);
  // partition slice starting at s+1
  p = s;
  for (i = s+1; i < s+n; i++)
    if (col[i] < x) {
      p++;
      swap(col, i, p);
	  swapd(val, i, p);
    }
  // swap pivot into place
  swap(col, s, p);
  swapd(val, s, p);
  // recurse into partition
  colQuickSort(col, val, s, p-s);
  colQuickSort(col, val, p+1, s+n-p-1);
}

void getRowIndex(mtxMatrix* A, int* d) {
    int curr_ind = 0;
	d[curr_ind] = 0;
	curr_ind++;
    for(int i = 0; i < A->NZ; i++) {
        if(A->Row[i] != A->Row[i+1]) {
            d[curr_ind] = i+1;
            curr_ind++;
        }
    }
}

int main(int argc, char *argv[])
{
    FILE *input_file, *output_file, *time_file;
    if (argc < 4)
	{
        fprintf(stderr, "Invalid input parameters\n");
		fprintf(stderr, "Usage: %s inputfile outputfile timefile \n", argv[0]);
		exit(1);
	}
    else    
    { 
        input_file = fopen(argv[1], "r");
        output_file = fopen(argv[2], "w");
        time_file = fopen(argv[3], "w");
        if(!input_file || !output_file || !time_file)
            exit(1);
    }

    MM_typecode matcode;
    if (mm_read_banner(input_file, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
  
    if (!mm_is_matrix(matcode) || !mm_is_real(matcode) ||
        !mm_is_coordinate(matcode) || !mm_is_symmetric(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    mtxMatrix inputMatrix;
    ReadMatrix(inputMatrix, input_file);

	quicksort(inputMatrix.Row, inputMatrix.Col, inputMatrix.Value, 0, inputMatrix.NZ);

    getRowIndex(&inputMatrix, inputMatrix.RowIndex);
	inputMatrix.RowIndex[inputMatrix.N] = inputMatrix.NZ;

	int *diag = new int[inputMatrix.N];
    
	for(int i = 0; i < inputMatrix.N + 1; i++) {
        for(int j = inputMatrix.RowIndex[i]; j < inputMatrix.RowIndex[i + 1]; j++) {
			colQuickSort(inputMatrix.Col + inputMatrix.RowIndex[i], inputMatrix.Value + inputMatrix.RowIndex[i], 0, inputMatrix.RowIndex[i + 1] - inputMatrix.RowIndex[i]);
		}
    }

	for(int i = 0; i < inputMatrix.N + 1; i++) {
        for(int j = inputMatrix.RowIndex[i]; j < inputMatrix.RowIndex[i + 1]; j++) {
			if (i == inputMatrix.Col[j]) diag[i] = j;
		}
    }

	

    for(int i = 0; i < inputMatrix.N + 1; i++) {
        printf("RowIndex[%i] = %i\n", i, inputMatrix.RowIndex[i]);
    }

     
    for(int i = 0; i < inputMatrix.N; i++) {
            printf("input[%i]= %lf\n", i, inputMatrix.Value[inputMatrix.RowIndex[i]]);
    }

	for(int i = 0; i < inputMatrix.N; i++) {
            printf("diag[%i]= %d\n", i, diag[i]);
    }

    ilu0(inputMatrix, inputMatrix.Value, diag);

	/*for(int i = 0; i < inputMatrix.NZ; i++) {
            printf("Value[%i]= %lf\n", i, inputMatrix.Value[i]);
    }*/
	WriteMatrix(inputMatrix, output_file, matcode);
   
    return 0;
}