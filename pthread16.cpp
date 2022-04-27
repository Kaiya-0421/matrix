#include <stdio.h>
#include <iostream>
//#include <Eigen/Dense>
#include <Eigen/Dense>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include <omp.h>
#include <pthread.h>
using Eigen::MatrixXd;
using namespace std;
using namespace Eigen;
#define M 1024
#define N 1024
//#define p_N 512 
#define bl 4 //ブロック数
void init(double Xmat[M][N]);                                                  //２次元配列の初期化関数
void *mymatrix(void *args);
void *mymatrix2(void *args); 
void *mymatrix3(void *args);
void *mymatrix4(void *args);
void *mymatrix5(void *args);
void *mymatrix6(void *args);
void *mymatrix7(void *args);
void *mymatrix8(void *args);
void *mymatrix9(void *args);
void *mymatrix10(void *args); 
void *mymatrix11(void *args);
void *mymatrix12(void *args);
void *mymatrix13(void *args);
void *mymatrix14(void *args);
void *mymatrix15(void *args);
void *mymatrix16(void *args);
        
void inimat(double m1[M][N], double m2[M][N], MatrixXd &mat1, MatrixXd &mat2); //自作行列関数用の配列同じランダム数を代入
// double sum(double result[M][N]);                                      //ずれの合計
double sum2(double result[M][N], MatrixXd &matresult2);
MatrixXd mat1 = MatrixXd::Random(M, N);
MatrixXd mat2 = MatrixXd::Random(M, N);
MatrixXd matresult;
double m1[M][N];
double m2[M][N];
double result[M][N];
int psum1=0,psum2=0,psum3=0,psum4=0,psum5=0,psum6=0,psum7=0,psum8=0;
int main(void)
{

  srand(1);
  double goukei;
  int i, j, k;
  struct timespec start, end;
  pthread_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11,t12, t13, t14, t15, t16;
  //配列初期化
  init(result);
  inimat(m1, m2, mat1, mat2);
  //時間計測開始
  clock_gettime(CLOCK_REALTIME, &start);
  //自作行列積関数
  pthread_create(&t1, NULL, mymatrix, NULL);
  pthread_create(&t2, NULL, mymatrix2, NULL);
  pthread_create(&t3, NULL, mymatrix3, NULL);
  pthread_create(&t4, NULL, mymatrix4, NULL);
  pthread_create(&t5, NULL, mymatrix5, NULL);
  pthread_create(&t6, NULL, mymatrix6, NULL);
  pthread_create(&t7, NULL, mymatrix7, NULL);
  pthread_create(&t8, NULL, mymatrix8, NULL);
  pthread_create(&t9, NULL, mymatrix9, NULL);
  pthread_create(&t10, NULL, mymatrix10, NULL);
  pthread_create(&t11, NULL, mymatrix11, NULL);
  pthread_create(&t12, NULL, mymatrix12, NULL);
  pthread_create(&t13, NULL, mymatrix13, NULL);
  pthread_create(&t14, NULL, mymatrix14, NULL);
  pthread_create(&t15, NULL, mymatrix15, NULL);
  pthread_create(&t16, NULL, mymatrix16, NULL);

  pthread_join(t1,NULL);
  pthread_join(t2,NULL);
  pthread_join(t3,NULL);
  pthread_join(t4,NULL);
  pthread_join(t5,NULL);
  pthread_join(t6,NULL);
  pthread_join(t7,NULL);
  pthread_join(t8,NULL);
  pthread_join(t9,NULL);
  pthread_join(t10,NULL);
  pthread_join(t11,NULL);
  pthread_join(t12,NULL);
  pthread_join(t13,NULL);
  pthread_join(t14,NULL);
  pthread_join(t15,NULL);
  pthread_join(t16,NULL);
  //行列積の終了時刻
  clock_gettime(CLOCK_REALTIME, &end);
  //外部ライブラリ行列積関数
  matresult = mat1 * mat2;
  // goukei = sum(result);
  goukei = sum2(result, matresult);
  printf("ずれの合計は%lf\n", goukei);
  //printf("平均のずれ率は%lf％\n", goukei * 100.0 / 9.0);

  
  double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
  printf("pthread時間は%lf\n", elasped_time);
  printf("aaa%daas",psum1+psum2+psum3+psum4+psum5+psum6+psum7+psum8);
  return 0;
}

//関数群
void init(double Xmat[M][N])
{
  int i, k;
  for (i = 0; i < M; i++)
  {
    for (k = 0; k < N - 4; k += 4)
    {
      Xmat[i][k] = 0.0;
      Xmat[i][k + 1] = 0.0;
      Xmat[i][k + 2] = 0.0;
      Xmat[i][k + 3] = 0.0;
    }
    for (; k < N; k++)
    {
      Xmat[i][k] = 0.0;
    }
  }
}
void *mymatrix(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 0; ii < 64; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];
                            ///
                           // psum1++;
                            //
                        }
                    }
                }
            }
        }
    }
    return NULL;
  
}
void *mymatrix2(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 64; ii < 128; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum2++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix3(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 128; ii < 192; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];
                            ///
                            //psum3++;
                            //
                        }
                    }
                }
            }
        }
    }
    return NULL;
  
}
void *mymatrix4(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 192; ii < 256; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum4++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix5(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 256; ii < 320; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum5++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix6(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 320; ii < 384; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum6++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix7(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 384; ii < 448; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                           //psum7++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix8(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 448; ii < 512; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum8++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix9(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 512; ii < 576; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];
                            ///
                           // psum1++;
                            //
                        }
                    }
                }
            }
        }
    }
    return NULL;
  
}
void *mymatrix10(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 576; ii < 640; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum2++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix11(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 640; ii < 704; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];
                            ///
                            //psum3++;
                            //
                        }
                    }
                }
            }
        }
    }
    return NULL;
  
}
void *mymatrix12(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 704; ii < 768; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum4++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix13(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 768; ii < 832; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum5++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix14(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 832; ii < 896; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                            //psum6++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix15(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 896; ii < 960; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                           //psum7++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}
void *mymatrix16(void *args)
{
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = 960; ii < 1024; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {

                for (i = ii; i < ii+bl; i++)
                {
                    for (j = jj; j < jj+bl; j++)
                    {
                        for (k = kk; k < kk+bl; k++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];

                            ///
                          //psum8++;
                            ///
                        }
                    }
                }
            }
        }
    }

  return NULL;

}


double sum2(double result[M][N], MatrixXd &matresult2)
{
  int i, k;
  double goukei = 0, goukei2 = 0, goukei3 = 0, goukei4 = 0;
  for (i = 0; i < M; i++)
  {
    for (k = 0; k < N - 4; k += 4)
    {
      goukei += fabs(result[i][k] - matresult2(i, k));
      goukei2 += fabs(result[i][k + 1] - matresult2(i, k + 1));
      goukei3 += fabs(result[i][k + 2] - matresult2(i, k + 2));
      goukei4 += fabs(result[i][k + 3] - matresult2(i, k + 3));
    }
    for (; k < N; k++)
    {
      goukei += fabs(result[i][k] - matresult2(i, k));
    }
  }
  goukei = goukei + goukei2 + goukei3 + goukei4;
  return goukei;
}
void inimat(double m1[M][N], double m2[M][N], MatrixXd &mat1, MatrixXd &mat2)
{
  int i, j;
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < N - 4; j += 4)
    {
      m1[i][j] = (mat1(i, j));
      m1[i][j + 1] = (mat1(i, j + 1));
      m1[i][j + 2] = (mat1(i, j + 2));
      m1[i][j + 3] = (mat1(i, j + 3));
      m2[i][j] = (mat2(i, j));
      m2[i][j + 1] = (mat2(i, j + 1));
      m2[i][j + 2] = (mat2(i, j + 2));
      m2[i][j + 3] = (mat2(i, j + 3));
    }
    //端数処理
    for (; j < N; j++)
    {
      m1[i][j] = (mat1(i, j));
      m2[i][j] = (mat2(i, j));
    }
  }
}
