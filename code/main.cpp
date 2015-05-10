#include <stdio.h>
#include "sparse.h"
#include "util.h"
#include "mmio.h"

const double EPSILON = 0.000001;

//int ilu0(mtxMatrix &A, double * luval, int * uptr)
//{
//     int j1, j2; // ������� ������� ������
//     int jrow; // ����� �������� �������
//     int k, j, jj; // �������� ������
//     int *iw = NULL; // ��������� ������
//     int jw;
//     double t1;
//
//     iw = new int[A.N];
//     memset(iw, 0, A.N * sizeof(int));
//
//     memcpy(luval, A.Value, A.NZ * sizeof(double));
//
//     for(k = 0; k < A.N; k++)
//     {
//        j1 = row[k];
//        j2 = row[k + 1];
//        for(j = j1; j < j2; j++)
//        {
//            iw[col[j]] = j;
//        }
//     }
//}

void getIndexesDiagonalElements(mtxMatrix* A, int* d) {
    if(A->Row[0] == A->Col[0])
        d[0] = 0;
    else
        d[0] = -1;

    int curr_ind = 1;
    for(int i = 1; i < A->NZ; i++) {
        if(A->Row[i] != A->Row[i+1]) {
            if (A->Row[i+1] == A->Col[i+1]) {
                d[curr_ind] = i+1;
            }
            else {
                d[curr_ind] = -1;
            }
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
    WriteMatrix(inputMatrix, output_file, matcode);

    int* d= new int[inputMatrix.N];
    getIndexesDiagonalElements(&inputMatrix, d);
    
    for(int i=0; i<inputMatrix.N; i++) {
        printf("d[%i] = %i\n", i, d[i]);
    }

     

    for(int i=0; i < inputMatrix.N; i++) {
        if (d[i] != -1)
            printf("input[%i]= %lf\n", i, inputMatrix.Value[d[i]]);
    }

    //ilu0(inputMatrix);
   
    return 0;
}