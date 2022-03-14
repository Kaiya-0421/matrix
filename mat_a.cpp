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


int main(void)
{
  constexpr int M=565;
  constexpr int N=565;
  srand(1);
  double goukei;
  int i, j, k;
  struct timespec start, end;
  MatrixXd mat1 = MatrixXd::Random(M, N);
  MatrixXd mat2 = MatrixXd::Random(M, N);
  MatrixXd matresult;
  //mat1/=10000000;
  //mat2/=10000000;
  clock_gettime(CLOCK_REALTIME,&start);
  double  m1[M][N];
  double  m2[M][N];
  double result[M][N];
  //自作関数行列積格納配列初期化
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < N; j++)
    {
      result[i][j] = 0.0;
    }
  }
  
  for (i = 0; i < M; i++)
  {
    for (j = 0; j < N; j++)
    {
      m1[i][j]=mat1.coeff(i,j);
      //std::cout << i << ", " << j << ", m1=" << m1[i][j] << ", mat1=" << mat1(i,j) << std::endl;
      m2[i][j]=(mat2(i,j));
    }
  }

//自作行列関数

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
  
 //外部ライブラリを使った行列積
  matresult = mat1 * mat2;
  //行列積の終了時刻
clock_gettime(CLOCK_REALTIME,&end);
  
  for (i = 0; i < M; i++)
  {
    for (k = 0; k < N; k++)
    {
      goukei += result[i][k]- matresult(i, k);
    }
  }
 // printf("result[24][24]=%lf\n\n",result[24][24]);
  printf("ずれの合計は%lf\n", goukei);
  printf("平均のずれ率は%lf％\n", goukei *100.0/ 9.0);
  // printf("%lf",std::fabs(result[17][17])-std::fabs(matresult(17,17)));
  /*

  std::cout << "Here is mat1*mat2:\n"
          << matresult(1,0)  <<std::endl;
          */
  
  //時間計測終了
  double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)*1e-9;
  printf("時間は%lfms\n",elasped_time);
  
}
