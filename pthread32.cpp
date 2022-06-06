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
#define M 2048
#define N 2048
#define p_N 32 //スレッド数
#define p_bl 64///p_bl*bl(ブロック数)=2048 //スレッドごとのループ数
#define bl 32 //ブロック数
void init(double Xmat[M][N]);                                                  //２次元配列の初期化関数
void mymatrix(int *argv);      //行列積関数 argvにスレッドごとのループ開始位置を渡す
void inimat(double m1[M][N], double m2[M][N], MatrixXd &mat1, MatrixXd &mat2); //自作行列関数用の配列同じランダム数を代入
// double sum(double result[M][N]);                                      //ずれの合計
double sum2(double result[M][N], MatrixXd &matresult2);
MatrixXd mat1 = MatrixXd::Random(M, N);
MatrixXd mat2 = MatrixXd::Random(M, N);
MatrixXd matresult;
double m1[M][N];
double m2[M][N];
double result[M][N];
int psum1=0;
int main(void)
{

  srand(1);
  double goukei;
  int i, j, k;
  struct timespec start, end;
  pthread_t t[p_N];
  //配列初期化
  init(result);
  inimat(m1, m2, mat1, mat2);
  //時間計測開始
  clock_gettime(CLOCK_REALTIME, &start);
  //自作行列積関数
  int p=0;
  for(p=0;p<p_N;p++){
    pthread_create(&t[p], NULL, (void *)mymatrix, (void *)&p);
  }
  for(p=0;p<p_N;p++){
      pthread_join(t[p],NULL);
  }




  //行列積の終了時刻
  clock_gettime(CLOCK_REALTIME, &end);
  //外部ライブラリ行列積関数
  matresult = mat1 * mat2;
  // goukei = sum(result);
  goukei = sum2(result, matresult);
  printf("ずれの合計は%lf\n", goukei);
  printf("平均のずれ率は%lf％\n", goukei * 100.0 / 9.0);

  
  double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
  printf("pthread時間は%lf\n", elasped_time);
  //printf("aaa%daas",psum1+psum2);
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
void mymatrix(int *argv)//argvにスレッドごとの引数,スレッド数32こなら0~31の数字が渡されるはず...
{ 
  printf("%d\n",argv);
  int p_roop=*argv;
 // printf("%d\n",*p_roop);
  int i, k, j;
  int ii, kk, jj;
  //#pragma omp parallel for num_threads(8)
      for (ii = p_roop*p_bl; ii < (p_roop+1)*p_bl; ii += bl)
    {
        for (jj = 0; jj < N; jj += bl)
        {
            for (kk = 0; kk < N; kk += bl)
            {
                for (i = ii; i < ii+bl; i++)
                {
                    for (k = kk; k < kk+bl; k++)
                    {
                        for (j = jj; j < jj+bl; j++)
                        {
                            result[i][j] += m1[i][k] * m2[k][j];
                           // psum1++;
                        }
                    }
                }
            }
      }
    }
  
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
