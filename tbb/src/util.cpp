#include "util.h"

void InitializeMatrix(int N, int NZ, mtxMatrix &mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ];
    mtx.Col = new int[NZ];
    mtx.Row = new int[NZ];
    mtx.RowIndex = new int[N + 1];
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

int ReadMatrix(mtxMatrix &mtx, FILE* inputfile)
{
    int M, N, nz;
    if (mm_read_mtx_crd_size(inputfile, &M, &N, &nz) !=0)
        exit(1);

    InitializeMatrix(N, nz, mtx);

    for (int i = 0; i < nz; i++)
    {
        fscanf(inputfile, "%i %i %lf", &(mtx.Row[i]), &(mtx.Col[i]), &(mtx.Value[i]));
        mtx.Col[i]--;
        mtx.Row[i]--;
    }
    fclose(inputfile);
    return 0;
}
