#include "flux.h"

int setICsubsup(int  nn,double * x,double * A,double * rho,double * V,double * T,double * P)
{
    int n;
    printf("Setting Initial Conditions for Problem Type: Subsonic-supersonic\n");
    //printf("           x/L   A/A*   rho/rho*   V/a0   T/T0   p/p0\n");

    for (n=0; n < nn; n++)
    {
      rho[n] = 1.0 - 0.3146*x[n];
      T[n]   = 1.0 - 0.2314*x[n];
      V[n]   = (0.1 + 1.09*x[n])*sqrt(T[n]);
      P[n]   = rho[n]*T[n];
      //printf("pt %3d: %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f\n",n,x[n],A[n],rho[n],V[n],T[n], P[n]);
    }
    return n;
}

int setICsubsonic(int  nn,double * x,double * A,double * rho,double * V,double * T,double * P)
{
    int n;
    printf("Setting Initial Conditions for Problem Type: Subsonic\n");
    //printf("           x/L   A/A*   rho/rho*   V/a0   T/T0   p/p0\n");

    for (n=0; n < nn; n++)
    {
      rho[n] = 1.0 - 0.023*x[n];
      T[n]   = 1.0 - 0.009333*x[n];
      V[n]   = 0.05 + 0.11*x[n];
      P[n]   = rho[n]*T[n];
      //printf("pt %3d: %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f\n",n,x[n],A[n],rho[n],V[n],T[n], P[n]);
    }
    return n;
}

