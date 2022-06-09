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
#define M 256
#define N 256
#define p_N 4 //スレッド数
#define p_bl 64///p_bl*p_N(スレッド数)=M //スレッドごとのループ間隔
#define bl 1 //ブロック数
void init(double Xmat[M][N]);                                                  //２次元配列の初期化関数
//void mymatrix(int *argv);      //行列積関数 argvにスレッドごとのループ開始位置を渡す
void *mymatrix(void *argv);      //行列積関数 argvにスレッドごとのループ開始位置を渡す
void inimat(double m1[M][N], double m2[M][N], MatrixXd &mat1, MatrixXd &mat2); //自作行列関数用の配列同じランダム数を代入
// double sum(double result[M][N]);                                      //ずれの合計
double sum(double result[M][N], MatrixXd &matresult);
MatrixXd mat1 = MatrixXd::Random(M, N);
MatrixXd mat2 = MatrixXd::Random(M, N);
MatrixXd matresult;
double m1[M][N];
double m2[M][N];
double result[M][N];
//int psum1=0;
int psum[p_N];
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
  int pNum[p_N];
  for(p=0;p<p_N;p++){
      pNum[p]=p;
    pthread_create(&t[p], NULL, mymatrix, (void *)&pNum[p]);
  }
  for(p=0;p<p_N;p++){
      pthread_join(t[p],NULL);
      psum[p_N]=0;
  }




  //行列積の終了時刻
  clock_gettime(CLOCK_REALTIME, &end);
  //外部ライブラリ行列積関数
  matresult = mat1 * mat2;

//////////////////////消せ




  goukei=0;
  // goukei = sum(result);
  goukei = sum(result, matresult);
  printf("ずれの合計は%lf\n", goukei);
  printf("平均のずれ率は%lf％\n", (double)goukei/(100.0*M));

  
  double elasped_time = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) * 1e-9;
  printf("pthread時間は%lf\n", elasped_time);
 // printf("aaa%d///%daas",psum[1],psum[21]);
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
void *mymatrix(void *argv)//argvにスレッドごとの引数,スレッド数32こなら0~31の数字が渡されるはず...
{
  //printf("%d\n",(int*) argv);
  //int p_roop=*argv;
  int p_roop = *((int*)argv);
  printf("スレッド依存数%d\n",p_roop);
  int i, k, j;
  int ii, kk, jj;

  //#pragma omp parallel for num_threads(8)


     for(i=p_roop;i<(p_roop+1)*p_bl;i++){
       for(j=0;j<M;j++){
         for(k=0;k<M;k++){
            result[i][j]+=m1[i][k]*m2[k][j];
            }
          //printf("[%d][%d]=%f\n",i,j,result[i][j]);
       }
     }



/*

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
                            printf("[%d][%d]=%f\n",i,j,result[i][j]);
                        }
                    }
                }
            }
         }
    }





    */
  return NULL;

}

double sum(double result[M][N], MatrixXd &matresult)
{
  int i, j;
  double goukei = 0.0;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4;j++)
    {
      printf("ex[%d][%d]=%f\n",i,j,matresult(i,j));
      printf("or[%d][%d]=%f\n",i,j,result[i,j]);
      goukei = result[i][j] - matresult(i,j);
      //printf("%lf",goukei);
    }
  }
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
