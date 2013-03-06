#include "flux.h"
#include "List.h"

/*
 *  Austin Flynt
 *  MacCormack 1D Euler Subsonic-Supersonic Nozzle
 *  Tue Mar  5 19:07:54 EST 2013
 */

int main(int argc, char * argv[])
{
  int i,j,k;
  int n,t;
  int nn;

  // File variables
  char filename[32];
  filename[0]='\0';
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
  double *rhon  = NULL;
  double *Tn    = NULL;
  double *Vn    = NULL;
  double *drdt  = NULL;
  double *dVdt  = NULL;
  double *dTdt  = NULL;
  double *drdtb = NULL;
  double *dVdtb = NULL;
  double *dTdtb = NULL;

  // Flow field variables
#if 0
  double Q1,Q2,Q3,Q4;
  double f1,f2,f3,f4;
  double fp1, fp2, fp3, fp4;
  double fm1, fm2, fm3, fm4;
#endif

  double nx,ny;
  double M  = 0.0;
  double C  = 0.5;
  double dt = 100.0;
  double dx = 0.1;
  double gma = 1.4;
  double alpha = 0.0;
  double mindt = 100.0;

#if 0
  double **DFP = NULL;
  double *DFPdata = NULL;
  double **DFM = NULL;
  double *DFMdata = NULL;

  alloc_2D(&DFP,&DFPdata, 4, 4);
  alloc_2D(&DFM,&DFMdata, 4, 4);
#endif

// ######################### Section Break ###########################
// feed nozzle shape



  // get Command line args
  parse_args(argc,argv,&nx,&ny,&M,&alpha,filename);

  // Calculate Q vector
  //initQ( &Q1, &Q2, &Q3, &Q4, M, &alpha);

// ######################### Section Break ###########################
// Read file
  /* Open 2D grid file for reading */
  if ((fp=fopen(filename,"r")) == NULL){
    printf("\nCould not open file <%s>\n",filename);
    exit(0);
  }

  /* Read number of grid points */
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of grid points
  sscanf(buff,"%d",&nn);
  printf("Number of grid points = %d\n",nn);

  if ((rho   = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for rho");exit(0);}
  if ((V     = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for V");exit(0);}
  if ((Mv    = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for Mv");exit(0);}
  if ((T     = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for T");exit(0);}
  if ((P     = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for P");exit(0);}
  if ((rhob  = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for rhob");exit(0);}
  if ((Tb    = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for Tb");exit(0);}
  if ((Vb    = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for Vb");exit(0);}
  if ((rhon  = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for rhon");exit(0);}
  if ((Tn    = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for Tn");exit(0);}
  if ((Vn    = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for Vn");exit(0);}
  if ((drdt  = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for drdt");exit(0);}
  if ((dVdt  = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for dVdt");exit(0);}
  if ((dTdt  = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for dTdt");exit(0);}
  if ((drdtb = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for drdtb");exit(0);}
  if ((dVdtb = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for dVdtb");exit(0);}
  if ((dTdtb = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for dTdtb");exit(0);}
  if ((x = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for X");
    exit(0);
  }
  if ((A = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for A");
    exit(0);
  }
  if ((lnA = (double*)malloc(nn*sizeof(double))) == NULL){
    printf("\nCould not allocate memory for lnA");
    exit(0);
  }



  // Read grid nodes x-locations
  for (n=0; n < nn; n++)
  {
    fgets(buff,bdim,fp);
    sscanf(buff,"%lg",&x[n]);
  }

  // Compute area
  for (n=0; n < nn; n++)
  {
    A[n] = 1.0 + 2.2*(x[n] - 1.5)*(x[n] - 1.5);
    lnA[n] = log(A[n]);
    //printf("pt %2d, ln(A) = %f\n", n, lnA[n]);
  }


  // Set IC
  //printf("Initial Conditions\n");
  //printf("           x/L   A/A*   rho/rho*   V/a0   T/T0\n");
  for (n=0; n < nn; n++)
  {
    rho[n] = 1.0 - 0.3146*x[n];
    T[n]   = 1.0 - 0.2314*x[n];
    V[n]   = (0.1 + 1.09*x[n])*sqrt(T[n]);
    P[n]   = rho[n]*T[n];
   // printf("pt %3d: %7.2f %7.3f %7.3f %7.3f %7.3f\n",n,x[n],A[n],rho[n],V[n],T[n]);
  }


  int maxiter = 1000;
  int ans = 0;

  // Solver Loop
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

    // BCs
    // Inflow
    V[0] = 2.0*V[1] - V[2];
    //printf("V1 = %f\n",V[0]);
    //printf("V2 = %f\n",V[1]);
    //printf("V3 = %f\n",V[2]);
    //printf("pt %3d: %7.5f, %7.5f, %7.5f, %7.5f\n",0,rho[0],V[0],T[0],P[0]);

    // Outflow
    rho[nn-1] = 2.0 * rho[nn-2] - rho[nn-3];
      V[nn-1] = 2.0 *   V[nn-2] -   V[nn-3];
      T[nn-1] = 2.0 *   T[nn-2] -   T[nn-3];
    //printf("pt %3d: %7.5f, %7.5f, %7.5f, %7.5f\n",i,rho[nn-1],V[nn-1],T[nn-1],P[nn-1]);


    //printf("Updated values\n");
    //printf("pt %3s: %7s, %7s, %7s, %7s, %7s\n","  i","    rho","      V","      T","     P","     M");

    //Update P,M
    for (i=0; i <= nn-1; i++)
    {
      P[i]   = rho[i]*T[i];
      Mv[i]  = V[i] / sqrt(T[i]);
    }


    printf("%4d, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f\n",n,rho[15],V[15],T[15],P[15],Mv[15]);

    // Continue?
    if(n== maxiter)
    {
      printf("Continue? : \n");
      scanf("%d",&ans);
      if(ans > 0)
      {
        maxiter += ans;
        ans = 0;
      }
    }
  }




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
  free(rhon ); rhon  = NULL;
  free(Tn   ); Tn    = NULL;
  free(Vn   ); Vn    = NULL;
  free(drdt ); drdt  = NULL;
  free(dVdt ); dVdt  = NULL;
  free(dTdt ); dTdt  = NULL;
  free(drdtb); drdtb = NULL;
  free(dVdtb); dVdtb = NULL;
  free(dTdtb); dTdtb = NULL;

  return 0;

}// end main
