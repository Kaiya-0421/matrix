#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cmath>
using Eigen::MatrixXd;
using namespace std;
using namespace Eigen;
#define M 1000
#define N 1000
void init(double Xmat[M][N]);                                         //２次元配列の初期化関数
void mymatrix(double m1[M][N], double m2[M][N], double result[M][N]); //行列積関数
MatrixXd mat1 = MatrixXd::Random(M, N);
MatrixXd mat2 = MatrixXd::Random(M, N);
MatrixXd matresult;
double m1[M][N];
double m2[M][N];
double result[M][N];
int main(void)
{
  // constexpr int M=565;
  // constexpr int N=565;
  srand(1);
  double goukei;
  int i, j, k;
  struct timespec start, end;


  //配列初期化
  init(result);
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < N; j++)
    {
      m1[i][j] = mat1.coeff(i, j);
      m2[i][j] = (mat2(i, j));
    }
  }
  //時間計測開始
  clock_gettime(CLOCK_REALTIME, &start);
  //自作行列積関数
  mymatrix(m1, m2, result);
  //行列積の終了時刻
  clock_gettime(CLOCK_REALTIME, &end);
  //外部ライブラリ行列積関数
  matresult = mat1 * mat2;


  for (i = 0; i < M; i++)
  {
    for (k = 0; k < N; k++)
    {
      goukei += fabs(result[i][k]) - fabs(matresult(i, k));
    }
  }

  printf("ずれの合計は%lf\n", goukei);
  printf("平均のずれ率は%lf％\n", goukei * 100.0 / 9.0);

  //時間計測終了
  double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
  printf("時間は%lfms\n", elasped_time);
  return 0;
}

//関数群
void init(double Xmat[M][N])
{
  int i, k;
  for (i = 0; i < M; i++)
  {
    for (k = 0; k < N; k++)
    {
      Xmat[i][k] = 0.0;
    }
  }
}
void mymatrix(double m1[M][N], double m2[M][N], double result[M][N])
{
  int i, k, j;
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < N; j++)
    {
      for (k = 0; k < N; k++)
      {
        result[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
}
