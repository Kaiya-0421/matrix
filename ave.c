#include <stdio.h>
#include <math.h>
int main()
{
    // double a[10]={};
    double a[10] = {};
    double ave = 0.0;
    int i;
    for (i = 0; i < 10; i++)
    {

        ave += a[i];
    }
    ave = ave / 10.0;
    printf("%lf", ave);
}