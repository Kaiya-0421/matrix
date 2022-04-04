#include <stdio.h>
#include <stdlib.h>
int main()
{
    double *matrix;
    int i, j, n, m;
    n = 100, m = 100;
    matrix = (double *) malloc(sizeof(double)* n *m);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j+=2)
        {
            matrix[i * m + j] = i * m + j;
            printf("%f\n", matrix[i * m + j]);
            matrix[i *m + (j+1)] = i * m + (j+1);
            printf("%f\n", matrix[i * m + (j+1)]);
        }
        /*for(;j<m;j++){
            matrix[i * m + j] = i * m + j;
            printf("%f\n", matrix[i * m + j]);
        }*/
    }
    free(matrix);
}