#include "flux.h"
#include "List.h"

/*
 *  Austin Flynt
 *  MacCormack 1D Euler Subsonic-Supersonic Nozzle
 *  Sat Mar 16 21:29:28 EDT 2013
 */

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

  // Conservative Variables
  double *U1    = NULL;
  double *U2    = NULL;
  double *U3    = NULL;

  double *Ub1   = NULL;
  double *Ub2   = NULL;
  double *Ub3   = NULL;

  double *dU1   = NULL;
  double *dU2   = NULL;
  double *dU3   = NULL;

  double *dUb1  = NULL;
  double *dUb2  = NULL;
  double *dUb3  = NULL;

  double *F1    = NULL;
  double *F2    = NULL;
  double *F3    = NULL;

  double *Fb1   = NULL;
  double *Fb2   = NULL;
  double *Fb3   = NULL;
  double *J2    = NULL;


  double nx,ny;
  double M  = 0.0;
  double alpha = 0.0;

  int maxiter = 1400;
  int ask = 0;
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

  // x and area
  safeAllocDouble(nn, &x     );
  safeAllocDouble(nn, &A     );
  safeAllocDouble(nn, &lnA   );

  // Primitives
  safeAllocDouble(nn, &rho   );
  safeAllocDouble(nn, &V     );
  safeAllocDouble(nn, &T     );

  // Dependent Primitives
  safeAllocDouble(nn, &P     );
  safeAllocDouble(nn, &Mv    );
  safeAllocDouble(nn, &Mv    );

  // predictor vars
  safeAllocDouble(nn, &rhob  );
  safeAllocDouble(nn, &Tb    );
  safeAllocDouble(nn, &Vb    );

  // partials
  safeAllocDouble(nn, &drdt  );
  safeAllocDouble(nn, &dVdt  );
  safeAllocDouble(nn, &dTdt  );

  // predictor partials
  safeAllocDouble(nn, &drdtb );
  safeAllocDouble(nn, &dVdtb );
  safeAllocDouble(nn, &dTdtb );

  // Allocate Conservative Vectors
  safeAllocDouble(nn, &U1    );
  safeAllocDouble(nn, &U2    );
  safeAllocDouble(nn, &U3    );

  safeAllocDouble(nn, &Ub1   );
  safeAllocDouble(nn, &Ub2   );
  safeAllocDouble(nn, &Ub3   );

  safeAllocDouble(nn, &dU1   );
  safeAllocDouble(nn, &dU2   );
  safeAllocDouble(nn, &dU3   );

  safeAllocDouble(nn, &dUb1  );
  safeAllocDouble(nn, &dUb2  );
  safeAllocDouble(nn, &dUb3  );

  safeAllocDouble(nn, &F1    );
  safeAllocDouble(nn, &F2    );
  safeAllocDouble(nn, &F3    );

  safeAllocDouble(nn, &Fb1   );
  safeAllocDouble(nn, &Fb2   );
  safeAllocDouble(nn, &Fb3   );

  safeAllocDouble(nn, &J2    );


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
  else if (ProblemType == 2)
    setICsubsonic(nn,x,A,rho,V,T,P);
  else
    setICcons(nn,x,A,rho,V,T,P,U1,U2,U3);

  exit(0);



  // ==========================================
  // Solver Loop
  // ==========================================
  q1dSolve(maxiter,nn,x,lnA,
           rho, V ,T ,
           rhob,Vb,Tb,
           drdt ,dVdt ,dTdt,
           drdtb,dVdtb,dTdtb,
           P, Mv, ProblemType, ask
           );

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

  free(U1  );    U1  = NULL;
  free(U2  );    U2  = NULL;
  free(U3  );    U3  = NULL;

  free(Ub1 );    Ub1 = NULL;
  free(Ub2 );    Ub2 = NULL;
  free(Ub3 );    Ub3 = NULL;

  free(dU1 );    dU1 = NULL;
  free(dU2 );    dU2 = NULL;
  free(dU3 );    dU3 = NULL;

  free(dUb1);    dUb1= NULL;
  free(dUb2);    dUb2= NULL;
  free(dUb3);    dUb3= NULL;

  free(F1  );    F1  = NULL;
  free(F2  );    F2  = NULL;
  free(F3  );    F3  = NULL;

  free(Fb1 );    Fb1 = NULL;
  free(Fb2 );    Fb2 = NULL;
  free(Fb3 );    Fb3 = NULL;

  free(J2  );    J2  = NULL;


  return 0;

}// end main


