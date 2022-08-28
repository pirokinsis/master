#include <iostream>

using namespace std;

extern "C"{void grad_(int n,double* b, double* d,double* s, double* q,double* wk,int* ix,double* dd);}
int main()
{

cout<<"jej";
int x=5;
double b[5][5];
double d[5][5];
double s[5][5];
double dd[5][5];
double q[5][5];
int ix[25];
int wk[25];
for (int i=0;i<5;i++)
    {
        for (int j=0;j<5;j++)
        {
        b[i][j]=1;
        d[i][j]=1;
        dd[i][j]=1;
        q[i][j]=1;
        ix[i*5+j]=1;
        wk[i*5+j]=0;
        } 
    }

grad_(x,(double *)b,(double *)d,(double *)s,(double *)q,(double *)wk,(int *)ix,(double *)dd);
}
