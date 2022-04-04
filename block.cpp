#include <stdio.h>
#include <time.h>
#define M  32 //ブロック幅
#define N 1024
double sum[M][M];
double m1[M][M];
double m2[M][M];
int main()
{

    int i, k, j;
    int ii, jj, kk;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    for (ii = 0; ii < N; ii += M)
    {
        for (jj = 0; jj < N; jj += M)
        {
            for (kk = 0; kk < N; kk += M)
            {

                for (i = 0; i < M; i++)
                {
                    for (j = 0; j < M; j++)
                    {
                        for (k = 0; k < M; k++)
                        {
                            sum[i][j] += m1[i][k] * m2[k][j];
                        }
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &end);
    double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("時間は%lfms\n", elasped_time);
}