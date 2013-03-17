#include "flux.h"
#include "List.h"

/*
 *  Austin Flynt
 *  MacCormack 1D Euler Subsonic-Supersonic Nozzle
 *  Sat Mar 16 21:29:28 EDT 2013
 */
int safeAllocDouble(const int nn, double ** V, char * name);

int main(int argc, char * argv[])
{
  int i,j,k;
  int n,t;
  int nn;

  // File variables
  char filename[32];
  filename[0]='\0';
  char outfile[32];
  outfile[0]='\0';
  FILE *fp = NULL;
  int bdim = 80;
  char buff[bdim];
  double *x = NULL;
  double *A = NULL;
  double *lnA = NULL;

  double *rho   = NULL;
  double *T     = NULL;
  double *V     = NULL;
  double *P     = NULL;
  double *Mv    = NULL;

  double *rhob  = NULL;
  double *Tb    = NULL;
  double *Vb    = NULL;

  double *drdt  = NULL;
  double *dVdt  = NULL;
  double *dTdt  = NULL;

  double *drdtb = NULL;
  double *dVdtb = NULL;
  double *dTdtb = NULL;


  double nx,ny;
  double M  = 0.0;
  double C  = 0.5;
  double dt = 100.0;
  double dx = 0.1;
  double gma = 1.4;
  double alpha = 0.0;
  double mindt = 100.0;
  double Pe = 0.93;

  int maxiter = 1400;
  int ask = 0;
  int ans = 0;
  int ProblemType = 0;

  // Flow field variables
#if 0
  double Q1,Q2,Q3,Q4;
  double f1,f2,f3,f4;
  double fp1, fp2, fp3, fp4;
  double fm1, fm2, fm3, fm4;
  double **DFP = NULL;
  double *DFPdata = NULL;
  double **DFM = NULL;
  double *DFMdata = NULL;

  alloc_2D(&DFP,&DFPdata, 4, 4);
  alloc_2D(&DFM,&DFMdata, 4, 4);
  // Calculate Q vector
  //initQ( &Q1, &Q2, &Q3, &Q4, M, &alpha);
#endif

// Parse Arguments
// ######################### Section Break ###########################


  // get Command line args
  parse_args(argc,argv,&maxiter,&ask,&ProblemType,&M,&alpha,filename,outfile);

  //printf("problem type = %d\n", ProblemType);
  //exit(0);


  // Open 2D grid file for reading
  if ((fp=fopen(filename,"r")) == NULL){
    printf("\nCould not open file <%s>\n",filename);
    exit(0);
  }

  // Read number of grid points
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of grid points
  sscanf(buff,"%d",&nn);
  printf("Number of grid points = %d\n",nn);

// Allocate memory
// ######################### Section Break ###########################

  safeAllocDouble(nn, &x     ,"x"    );
  safeAllocDouble(nn, &A     ,"A"    );
  safeAllocDouble(nn, &lnA   ,"lnA"  );
  safeAllocDouble(nn, &rho   ,"rho"  );
  safeAllocDouble(nn, &V     ,"V"    );
  safeAllocDouble(nn, &Mv    ,"Mv"   );
  safeAllocDouble(nn, &Mv    ,"Mv"   );
  safeAllocDouble(nn, &T     ,"T"    );
  safeAllocDouble(nn, &P     ,"P"    );
  safeAllocDouble(nn, &rhob  ,"rhob" );
  safeAllocDouble(nn, &Tb    ,"Tb"   );
  safeAllocDouble(nn, &Vb    ,"Vb"   );
  safeAllocDouble(nn, &drdt  ,"drdt" );
  safeAllocDouble(nn, &dVdt  ,"dVdt" );
  safeAllocDouble(nn, &dTdt  ,"dTdt" );
  safeAllocDouble(nn, &drdtb ,"drdtb");
  safeAllocDouble(nn, &dVdtb ,"dVdtb");
  safeAllocDouble(nn, &dTdtb ,"dTdtb");


  // Read grid nodes x-locations and Area
  for (n=0; n < nn; n++)
  {
    fgets(buff,bdim,fp);
    sscanf(buff,"%lg %lg",&x[n],&A[n]);
    sscanf(buff,"%lg",&x[n]);
    //printf("From file: x = %f Area = %f\n", x[n], A[n]);
  }
  fclose(fp); fp = NULL;


  // Compute Log area
  for (n=0; n < nn; n++)
  {
    lnA[n] = log(A[n]);
    //printf("pt %2d, ln(A) = %f\n", n, lnA[n]);
    //A[n] = 1.0 + 2.2*(x[n] - 1.5)*(x[n] - 1.5);
    //printf("%.15f %.15f\n", x[n], A[n]);
  }


  // Set IC
  if(ProblemType == 1)
    setICsubsup(nn,x,A,rho,V,T,P);
  else
    setICsubsonic(nn,x,A,rho,V,T,P);



  // ==========================================
  // Solver Loop
  // ==========================================
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
      printf("How many more iterations?:>  ");
      scanf("%d",&ans);
      if(ans > 0)
      {
        maxiter += ans;
        ans = 0;
      }
    }
  }
  // ======================
  //  end solver iteration
  // ======================

  // Print soln to screen
  printSoln(nn,x,A,rho,V,T,P,Mv);



  // Open file for write
  if ((fp=fopen(outfile,"w")) == NULL){
    printf("\nCould not open file <%s>\n",outfile);
    exit(0);
  }

  // Write solution to file
  printf("writing to file: %s\n",outfile);
  writeSoln(fp,nn,x,A,rho,V,T,P,Mv);
  fclose(fp); fp = NULL;


  // Open file for write
  if ((fp=fopen("p.dat","w")) == NULL){printf("\nCould not open file <%s>\n","p.dat");exit(0);}
  write1var(fp, nn,  x,  P);
  fclose(fp); fp = NULL;

  if ((fp=fopen("r.dat","w")) == NULL){printf("\nCould not open file <%s>\n","r.dat");exit(0);}
  write1var(fp, nn,  x,  rho);
  fclose(fp); fp = NULL;

  if ((fp=fopen("T.dat","w")) == NULL){printf("\nCould not open file <%s>\n","T.dat");exit(0);}
  write1var(fp, nn,  x,  T);
  fclose(fp); fp = NULL;

  if ((fp=fopen("v.dat","w")) == NULL){printf("\nCould not open file <%s>\n","v.dat");exit(0);}
  write1var(fp, nn,  x,  V);
  fclose(fp); fp = NULL;


  // Free memory
  free(x);         x = NULL;
  free(A);         A = NULL;
  free(lnA);     lnA = NULL;
  free(rho);     rho = NULL;
  free(T);         T = NULL;
  free(V);         V = NULL;
  free(P);         P = NULL;
  free(Mv);       Mv = NULL;
  free(rhob);   rhob = NULL;
  free(Tb);       Tb = NULL;
  free(Vb   ); Vb    = NULL;
  free(drdt ); drdt  = NULL;
  free(dVdt ); dVdt  = NULL;
  free(dTdt ); dTdt  = NULL;
  free(drdtb); drdtb = NULL;
  free(dVdtb); dVdtb = NULL;
  free(dTdtb); dTdtb = NULL;

  return 0;

}// end main



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

int safeAllocDouble(const int nn, double ** V, char * name)
{
  if (((*V)   = (double*)malloc(nn*sizeof(double))) == NULL){
    //printf("\nCould not allocate memory for %s", varname);
    printf("\nCould not allocate memory for %s", name);
    exit(0);
  }
}
