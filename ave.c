#include <stdio.h>
#include <math.h>
int main()
{
    // double a[10]={};
    double a[10] = {2.892361,2.888346,2.861979,2.872691,2.885919,2.886732,2.871982,2.877001,3.027935,2.885614};
    double ave = 0.0;
    int i;
    for (i = 0; i < 10; i++)
    {

        ave += a[i];
    }
    ave = ave / 10.0;
    printf("%lf", ave);
}