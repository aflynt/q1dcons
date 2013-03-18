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

// Set Conservative IC
int setICcons(int  nn,double * x,double * A,double * rho,double * V,double * T,double * P,
              double * U1, double * U2, double * U3)
{
    int n;
    double mdot = 0.59;

    printf(" In CONS IC\n");
    printf("Setting Initial Conditions for Problem Type: Conservative sub-super\n");
    printf("           x/L   A/A*   rho/rho*   V/a0   T/T0    U1      U2     U3\n");

    // Loop over nodes
    for (n=0; n < nn; n++)
    {
      if((x[n] -0.5) < 1.e-10) // x less than 0.5
      {
        rho[n] = 1.0;
          T[n] = 1.0;
      }

      else if((x[n] - 1.5) < 1.e-10) // x less than 1.5
      {
        rho[n] = 1.0 - 0.366*(x[n] - 0.5);
        T[n]   = 1.0 - 0.167*(x[n] - 0.5);
      }

      else // x greater than 1.5
      {
        rho[n] = 0.634 - 0.3879*(x[n] - 1.5);
        T[n]   = 0.833 - 0.3507*(x[n] - 1.5);
      }

      V[n]   = mdot/(rho[n]*A[n]);
      P[n]   = rho[n]*T[n];

      // U vars
      U1[n] = rho[n]*A[n];
      U2[n] = rho[n]*A[n]*V[n];
      U3[n] = rho[n]* (T[n]/gm1 + g/2.0*V[n]*V[n])*A[n];


      printf("pt %3d: %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
                n,x[n],A[n],rho[n],V[n],T[n], U1[n], U2[n], U3[n]);
    }
    return n;
}

// QQ
int  q1dSolve(int maxiter,const int nn, double * x, double * lnA,
            double * rho,   double * V ,    double * T ,
            double * rhob,  double * Vb,    double * Tb,
            double * drdt , double * dVdt , double * dTdt,
            double * drdtb, double * dVdtb, double * dTdtb,
            double * P,     double * Mv, const int ProblemType, int ask)
{
  int i, n;
  int ans = 0;
  double gma = 1.4;
  double mindt = 100.0;
  double dt = 100.0;
  double dx = 0.1;
  double C  = 0.5;
  double Pe = 0.93;

  for (n = 1; n <= maxiter; n++)
  {


    // Get Predictor Partials
    //printf("Predictor Partials\n");
    for (i=0; i < nn-1; i++)
    {
      drdt[i]  = -rho[i] * (V[i+1] - V[i])/dx
                - rho[i] * V[i] * (lnA[i+1] - lnA[i])/dx
                -   V[i] * (rho[i+1] - rho[i])/dx;

      dVdt[i]  = -V[i] * (V[i+1] - V[i])/dx
                 - 1.0/gma*((T[i+1] - T[i])/dx + T[i]/rho[i]*(rho[i+1] - rho[i])/dx);

      dTdt[i]  = -V[i] * (T[i+1] - T[i])/dx
                 - (gma-1.0)*T[i]*((V[i+1] - V[i])/dx + V[i]*(lnA[i+1] - lnA[i])/dx);
      //printf("pt %3d: %7.5f %7.5f %7.5f\n",i,drdt[i],dVdt[i],dTdt[i]);
    }


    // Get min timestep
    for (i=1; i < nn-1; i++)
    {
      dt = C*dx/(sqrt(T[i]) + V[i]);
      mindt = min(mindt,dt);
    }
    //printf("minimum dt = %f\n", mindt);



    // Calc Predictor values
    //printf("Calculating Predictor values\n");
    for (i=0; i < nn-1; i++)
    {
      rhob[i] = rho[i] + drdt[i]*mindt;
        Vb[i] =   V[i] + dVdt[i]*mindt;
        Tb[i] =   T[i] + dTdt[i]*mindt;
      //printf("pt %3d: %7.5f %7.5f %7.5f\n",i,rhob[i],Vb[i],Tb[i]);
    }

    //CORRECTOR STEP
    //
    //
    //Calc Corrector Partials
    //printf("Corrector Partials\n");
    //printf("pt %3s: %7s, %7s, %7s\n","  i","   drho","     dV","     dT");
    for (i=1; i < nn-1; i++)
    {
      drdtb[i] =  -rhob[i]*(Vb[i] - Vb[i-1])/dx
                 - rhob[i]*Vb[i]*(lnA[i]-lnA[i-1])/dx
                 -Vb[i]*(rhob[i] - rhob[i-1])/dx;
      dVdtb[i] = -Vb[i]*(Vb[i] - Vb[i-1])/dx
                 -1.0/gma *
                  ((Tb[i] - Tb[i-1])/dx  + Tb[i]/rhob[i]*(rhob[i]-rhob[i-1])/dx);
      dTdtb[i] = -Vb[i]*(Tb[i] - Tb[i-1])/dx
                 - (gma-1.0)*Tb[i] *
                   ((Vb[i] - Vb[i-1])/dx + Vb[i]*(lnA[i]-lnA[i-1])/dx);
      //printf("pt %3d: %7.5f, %7.5f, %7.5f\n",i,drdtb[i],dVdtb[i],dTdtb[i]);
    }

    // Average partials
    //printf("Average Partials\n");
    //printf("pt %3s: %7s, %7s, %7s\n","  i","av drho","av   dV","av   dT");
    for (i=1; i < nn-1; i++)
    {
      drdt[i] = 0.5*(drdt[i] + drdtb[i]);
      dVdt[i] = 0.5*(dVdt[i] + dVdtb[i]);
      dTdt[i] = 0.5*(dTdt[i] + dTdtb[i]);
      //printf("pt %3d: %7.5f, %7.5f, %7.5f\n",i,drdt[i],dVdt[i],dTdt[i]);
    }

    // Calc Corrector values
    for (i=1; i < nn-1; i++)
    {
      rho[i] =  rho[i] + drdt[i]*mindt;
      V[i]   =    V[i] + dVdt[i]*mindt;
      T[i]   =    T[i] + dTdt[i]*mindt;
    }

    // CORRECTOR COMPLETE

    // =========
    // BCs
    // =========
    // Inflow
    // rho, T remain constant
    V[0] = 2.0*V[1] - V[2];

    // Outflow
    if(ProblemType == 1) // subsonic - supersonic
    {
    rho[nn-1] = 2.0 * rho[nn-2] - rho[nn-3]; //extrapolate rho'
      V[nn-1] = 2.0 *   V[nn-2] -   V[nn-3]; //extrapolate V'
      T[nn-1] = 2.0 *   T[nn-2] -   T[nn-3]; //extrapolate T'
    }
    else // subsonic P = specified for subsonic
    {
      P[nn-1] = Pe;                          //specify outlet pressure
      T[nn-1] = 2.0 *   T[nn-2] -   T[nn-3]; //extrapolate T
    rho[nn-1] = P[nn-1] / T[nn-1];           //extrapolate rho'= P/T
      V[nn-1] = 2.0 *   V[nn-2] -   V[nn-3]; //extrapolate V
    }
    //printf("pt %3d: %7.5f, %7.5f, %7.5f, %7.5f\n",i,rho[nn-1],V[nn-1],T[nn-1],P[nn-1]);


    //printf("Updated values\n");
    //printf("pt %3s: %7s, %7s, %7s, %7s, %7s\n","  i","    rho","      V","      T","     P","     M");

    //Update P,M
    if(ProblemType == 1) // subsonic - supersonic
    {
      for (i=0; i <= nn-1; i++)
      {
        P[i]   = rho[i]*T[i];       // Compute dependent P
        Mv[i]  = V[i] / sqrt(T[i]); // Compute dependent M
      }
    }
    else // subsonic
    {
      for (i=0; i < nn-1; i++)
        P[i]   = rho[i]*T[i];       // only update internal Pressure values

      for (i=0; i <= nn-1; i++)
        Mv[i]  = V[i] / sqrt(T[i]); // compute M across whole field
    }


    // Print values at nozzle throat
    printf(" %4d", n);
    if (n %20 == 0)
      putchar('\n');
    //printf("%4d, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n",
    //          n,rho[15],V[15],T[15],P[15],Mv[15]);

    // Continue?
    if(n == maxiter && ask)
    {
      printf("\nHow many more iterations?:>  ");
      scanf("%d",&ans);
      if(ans > 0)
      {
        maxiter += ans;
        ans = 0;
      }
    }
  }

  return maxiter;
} // end of solver


// XX
int  q1dSolveCons(int maxiter,const int nn, double * x, double * A, double * lnA,
            double * rho,   double * V ,    double * T ,
            double * rhob,  double * Vb,    double * Tb,
            double * drdt,  double * dVdt,  double * dTdt,
            double * drdtb, double * dVdtb, double * dTdtb,
            double * U1,    double * U2,    double * U3,
            double * Ub1,   double * Ub2,   double * Ub3,
            double * dU1,   double * dU2,   double * dU3,
            double * dUb1,  double * dUb2,  double * dUb3,
            double * F1,    double * F2,    double * F3,
            double * Fb1,   double * Fb2,   double * Fb3,
            double * J2,    double * P,     double * Mv,
            const int ProblemType, int ask)

{
  int i, n;
  int ans = 0;
  double gma = 1.4;
  double mindt = 100.0;
  double dt = 100.0;
  double dx = 0.1;
  double C  = 0.5;
  double Pe = 0.93;

  for (n = 1; n <= maxiter; n++)
  {

    //Calc Fluxes
    printf("\nFluxes =\n");
    printf("pt %3s: %7s %7s %7s %7s\n","itr", "    F1", "    F2", "    F3", "    J2");
    for (i=0; i < nn-1; i++)
    {
      F1[i] = U2[i];
      F2[i] = U2[i]*U2[i]/U1[i] + gm1/g*(U3[i] - g/2.0*U2[i]*U2[i]/U1[i]);
      F3[i] = g*U2[i]*U3[i]/U1[i] - g*gm1/2.0*U2[i]*U2[i]*U2[i]/(U1[i]*U1[i]);
      J2[i] = 1/g*rho[i]*T[i]*(A[i+1] - A[i])/dx;
      printf("pt %3d: %7.2f %7.3f %7.3f %7.3f\n",i, F1[i], F2[i], F3[i], J2[i]);
    }
    return maxiter;



    // Get Predictor Partials
    //printf("Predictor Partials\n");
    for (i=0; i < nn-1; i++)
    {
      drdt[i]  = -rho[i] * (V[i+1] - V[i])/dx
                - rho[i] * V[i] * (lnA[i+1] - lnA[i])/dx
                -   V[i] * (rho[i+1] - rho[i])/dx;

      dVdt[i]  = -V[i] * (V[i+1] - V[i])/dx
                 - 1.0/gma*((T[i+1] - T[i])/dx + T[i]/rho[i]*(rho[i+1] - rho[i])/dx);

      dTdt[i]  = -V[i] * (T[i+1] - T[i])/dx
                 - (gma-1.0)*T[i]*((V[i+1] - V[i])/dx + V[i]*(lnA[i+1] - lnA[i])/dx);
      //printf("pt %3d: %7.5f %7.5f %7.5f\n",i,drdt[i],dVdt[i],dTdt[i]);
    }


    // Get min timestep
    for (i=1; i < nn-1; i++)
    {
      dt = C*dx/(sqrt(T[i]) + V[i]);
      mindt = min(mindt,dt);
    }
    //printf("minimum dt = %f\n", mindt);



    // Calc Predictor values
    //printf("Calculating Predictor values\n");
    for (i=0; i < nn-1; i++)
    {
      rhob[i] = rho[i] + drdt[i]*mindt;
        Vb[i] =   V[i] + dVdt[i]*mindt;
        Tb[i] =   T[i] + dTdt[i]*mindt;
      //printf("pt %3d: %7.5f %7.5f %7.5f\n",i,rhob[i],Vb[i],Tb[i]);
    }

    //CORRECTOR STEP
    //
    //
    //Calc Corrector Partials
    //printf("Corrector Partials\n");
    //printf("pt %3s: %7s, %7s, %7s\n","  i","   drho","     dV","     dT");
    for (i=1; i < nn-1; i++)
    {
      drdtb[i] =  -rhob[i]*(Vb[i] - Vb[i-1])/dx
                 - rhob[i]*Vb[i]*(lnA[i]-lnA[i-1])/dx
                 -Vb[i]*(rhob[i] - rhob[i-1])/dx;
      dVdtb[i] = -Vb[i]*(Vb[i] - Vb[i-1])/dx
                 -1.0/gma *
                  ((Tb[i] - Tb[i-1])/dx  + Tb[i]/rhob[i]*(rhob[i]-rhob[i-1])/dx);
      dTdtb[i] = -Vb[i]*(Tb[i] - Tb[i-1])/dx
                 - (gma-1.0)*Tb[i] *
                   ((Vb[i] - Vb[i-1])/dx + Vb[i]*(lnA[i]-lnA[i-1])/dx);
      //printf("pt %3d: %7.5f, %7.5f, %7.5f\n",i,drdtb[i],dVdtb[i],dTdtb[i]);
    }

    // Average partials
    //printf("Average Partials\n");
    //printf("pt %3s: %7s, %7s, %7s\n","  i","av drho","av   dV","av   dT");
    for (i=1; i < nn-1; i++)
    {
      drdt[i] = 0.5*(drdt[i] + drdtb[i]);
      dVdt[i] = 0.5*(dVdt[i] + dVdtb[i]);
      dTdt[i] = 0.5*(dTdt[i] + dTdtb[i]);
      //printf("pt %3d: %7.5f, %7.5f, %7.5f\n",i,drdt[i],dVdt[i],dTdt[i]);
    }

    // Calc Corrector values
    for (i=1; i < nn-1; i++)
    {
      rho[i] =  rho[i] + drdt[i]*mindt;
      V[i]   =    V[i] + dVdt[i]*mindt;
      T[i]   =    T[i] + dTdt[i]*mindt;
    }

    // CORRECTOR COMPLETE

    // =========
    // BCs
    // =========
    // Inflow
    // rho, T remain constant
    V[0] = 2.0*V[1] - V[2];

    // Outflow
    if(ProblemType == 1) // subsonic - supersonic
    {
    rho[nn-1] = 2.0 * rho[nn-2] - rho[nn-3]; //extrapolate rho'
      V[nn-1] = 2.0 *   V[nn-2] -   V[nn-3]; //extrapolate V'
      T[nn-1] = 2.0 *   T[nn-2] -   T[nn-3]; //extrapolate T'
    }
    else // subsonic P = specified for subsonic
    {
      P[nn-1] = Pe;                          //specify outlet pressure
      T[nn-1] = 2.0 *   T[nn-2] -   T[nn-3]; //extrapolate T
    rho[nn-1] = P[nn-1] / T[nn-1];           //extrapolate rho'= P/T
      V[nn-1] = 2.0 *   V[nn-2] -   V[nn-3]; //extrapolate V
    }
    //printf("pt %3d: %7.5f, %7.5f, %7.5f, %7.5f\n",i,rho[nn-1],V[nn-1],T[nn-1],P[nn-1]);


    //printf("Updated values\n");
    //printf("pt %3s: %7s, %7s, %7s, %7s, %7s\n","  i","    rho","      V","      T","     P","     M");

    //Update P,M
    if(ProblemType == 1) // subsonic - supersonic
    {
      for (i=0; i <= nn-1; i++)
      {
        P[i]   = rho[i]*T[i];       // Compute dependent P
        Mv[i]  = V[i] / sqrt(T[i]); // Compute dependent M
      }
    }
    else // subsonic
    {
      for (i=0; i < nn-1; i++)
        P[i]   = rho[i]*T[i];       // only update internal Pressure values

      for (i=0; i <= nn-1; i++)
        Mv[i]  = V[i] / sqrt(T[i]); // compute M across whole field
    }


    // Print values at nozzle throat
    printf(" %4d", n);
    if (n %20 == 0)
      putchar('\n');
    //printf("%4d, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n",
    //          n,rho[15],V[15],T[15],P[15],Mv[15]);

    // Continue?
    if(n == maxiter && ask)
    {
      printf("\nHow many more iterations?:>  ");
      scanf("%d",&ans);
      if(ans > 0)
      {
        maxiter += ans;
        ans = 0;
      }
    }
  }

  return maxiter;
} // end of solver
