#include "util.h"
#include "mmio.h"
#include "timer.h"
#include <fstream>
#include "memory.h"

const double EPSILON = 0.000001;


int ilu0(mtxMatrix &A, double * luval, int * uptr) {
    int jrow = 0;
	int jstatus = 0;
    int *iw = NULL;
    int jw = 0;
	int jj = 0;
    double t1;

	iw = new int[A.N];
	memset(iw, 0, A.N * sizeof(int));
    memcpy(luval, A.Value, A.RowIndex[A.N] * sizeof(double));

    bool flag1 = true;
		
	for (int k = 0; (k < A.N) && flag1; k++) {
		int j1 = A.RowIndex[k];
		int j2 = A.RowIndex[k + 1];


		//if (!flag1) continue;
		printf("k= %d, j1= %d, j2= %d \n", k, j1 , j2);
		
		for (int j = j1; j < j2; j++) {
			iw[A.Col[j]] = j;
		}

		int j = j1;
		int jstatus = j1;
		bool flag2;

		#pragma omp parallel private(j, jj, jrow, t1, jw, flag2)
		{		
			flag2 = (j1 < j2) && (A.Col[j1] < k);
			#pragma omp for
			for (j = j1; j < j2; j++) {
				flag2 = (A.Col[j] < k);
				if (!flag2) continue;
				printf("k=%d , j= %d \n", k, j);
			
				jrow = A.Col[j];
				t1 = luval[j] / luval[uptr[jrow]];
				luval[j] = t1;

				for (jj = uptr[jrow] + 1; jj < A.RowIndex[jrow + 1]; jj++) {
					jw = iw[A.Col[jj]];
					if (jw != 0) {
						luval[jw] = luval[jw] - t1 * luval[jj];
					}
				}
				flag2 = (j + 1 < j2) && (A.Col[j + 1] < k);
				if (!flag2) jstatus = j + 1;
			}
		}
		j = jstatus;
		jrow = A.Col[j];
		uptr[k] = j;

		flag1 = !((jrow != k) || (fabs(luval[j]) < EPSILON)) && (k + 1 < A.N);

		for (j = j1; j < j2; j++) {
			iw[A.Col[j]] = 0;
		}
	}

    delete[] iw;

    return 0;
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

int main(int argc, char *argv[]) {
    FILE *input_file, *output_file, *time_file;
    if (argc < 4) {
        fprintf(stderr, "Invalid input parameters\n");
        fprintf(stderr, "Usage: %s inputfile outputfile timefile \n", argv[0]);
        exit(1);
    }
    else {
        input_file = fopen(argv[1], "r");
        output_file = fopen(argv[2], "w");
        time_file = fopen(argv[3], "w");
        if (!input_file || !output_file || !time_file)
            exit(1);
    }

    MM_typecode matcode;
    if (mm_read_banner(input_file, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (!mm_is_matrix(matcode) || !mm_is_real(matcode) || !mm_is_coordinate(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    mtxMatrix inputMatrix, fullMatrix;
    ReadMatrix(inputMatrix, input_file);
    Timer timer;

    getRowIndex(&inputMatrix, inputMatrix.RowIndex);
    inputMatrix.RowIndex[inputMatrix.N] = inputMatrix.NZ;

    if (mm_is_symmetric(matcode)) {
        InitializeMatrix(inputMatrix.N, 2 * inputMatrix.NZ - inputMatrix.N, fullMatrix);
        TriangleToFull(&inputMatrix, &fullMatrix);
        FreeMatrix(inputMatrix);
    }
    else {
        fullMatrix = inputMatrix;
    }

    int *diag = new int[fullMatrix.N];

    for (int i = 0; i < fullMatrix.N + 1; i++) {
        for (int j = fullMatrix.RowIndex[i]; j < fullMatrix.RowIndex[i + 1]; j++) {
            if (i == fullMatrix.Col[j]) diag[i] = j;
        }
    }

//    for (int i = 0; i < fullMatrix.N + 1; i++) {
//        printf("RowIndex[%i] = %i\n", i, fullMatrix.RowIndex[i]);
//    }
//
//
//    for (int i = 0; i < fullMatrix.N; i++) {
//        printf("input[%i]= %lf\n", i, fullMatrix.Value[inputMatrix.RowIndex[i]]);
//    }
//
//    for (int i = 0; i < fullMatrix.N; i++) {
//        printf("diag[%i]= %d\n", i, diag[i]);
//    }

    timer.start();
    ilu0(fullMatrix, fullMatrix.Value, diag);
    timer.stop();

    std::ofstream timeLog;
    timeLog.open(argv[3]);
    timeLog << timer.getElapsed();

    WriteFullMatrix(fullMatrix, output_file, matcode);

    FreeMatrix(fullMatrix);

    return 0;
}
