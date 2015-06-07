class WawePropagation
{
    int type;//0-increase,1-stable;2-decrease
    double **M, **B, *dm, x_step, y_step, omega;
    int wawe, k, nx, ny;

public:
    WawePropagation(double **tM, double **tB, int twawe, double *tdm, int ttype, \
        int tk, double tx_step, double ty_step, int tnx, int tny, double tomega) :
                    M(tM), B(tB), wawe(twawe), dm(tdm), type(ttype), k(tk), x_step(tx_step), y_step(ty_step), nx(tnx), ny(tny), omega(tomega) {}
    void operator()( const blocked_range<int>& r ) const
    {
        int begin=r.begin(), end=r.end(),j;
        double temp, d;
        double x_st_2=x_step*x_step, y_st_2 = y_step*y_step;
        if (type==0)
            for ( int i=begin; i<end; i++ ) // i 1:wawe
            {
                j = wawe + 1 - i; // j wawe:1
                temp = M[i][j];
                M[i][j] =  ( 1 - omega)*M[i][j] + omega * 0.5 * (	y_st_2*(M[i-1][j] +
                                    M[i+1][j]) +
                                    x_st_2*(M[i][j-1] +
                                    M[i][j+1]) -
                                    x_st_2*y_st_2* B[i][j]
                                    ) / (x_st_2 + y_st_2);
                //printf("M[%i][%i]=%f\n",i,j,M[i][j]);
                d = fabs(temp-M[i][j]);
                if ( dm[wawe] < d ) dm[wawe] = d;
            }
        else if (type == 1)
            for ( int i=begin; i<end; i++  ) // i 1:wawe
            {
                j = wawe + k + 1 - i; // j wawe+k:k+1
                temp = M[i][j];
                M[i][j] =  ( 1 - omega)*M[i][j] + omega * 0.5 * (	y_st_2*(M[i-1][j] +
                                    M[i+1][j]) +
                                    x_st_2*(M[i][j-1] +
                                    M[i][j+1]) -
                                    x_st_2*y_st_2* B[i][j]
                                    ) / (x_st_2 + y_st_2);
                //printf("M[%i][%i]=%f\n",i,j,M[i][j]);
                d = fabs(temp-M[i][j]);
                if ( dm[wawe] < d ) dm[wawe] = d;
            }
        else if (type == 2)
            for ( int i=begin; i<end; i++  ) // i nx-1-wawe:nx-2
            {
                j = ny+nx-wawe-3-i; //  j ny-2:ny-1-wawe
                temp = M[i][j];
                M[i][j] =  ( 1 - omega)*M[i][j] + omega * 0.5 * (	y_st_2*(M[i-1][j] +
                                    M[i+1][j]) +
                                    x_st_2*(M[i][j-1] +
                                    M[i][j+1]) -
                                    x_st_2*y_st_2* B[i][j]
                                    ) / (x_st_2 + y_st_2);
                //printf("M[%i][%i]=%f\n",i,j,M[i][j]);
                d = fabs(temp-M[i][j]);
                if ( dm[wawe] < d ) dm[wawe] = d;
            }
    }
};
