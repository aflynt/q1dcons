#include "flux.h"

int help(void) {
  printf("Usage: ./flux  [-nxy <nx> <ny>] [-Q <Q1 Q2 Q3 Q4>] [-omega <omega>] [-P <filename>]\n");
  printf("\t-nxy: follow by area components in x-dir and y-dir\n");
  printf("\t-Q:   follow by Q vector components\n");
  printf("\t-omega: SOR relaxation factor (0 < omega < 2.0)\n");
  printf("\t-P: Print T to file <filename>\n");

  return 1;
}

// Open FILE 'filename' pointed to by 'fp'
int open_file(char * filename, FILE ** fp)
{
  if (((*fp)=fopen(filename,"r")) == NULL)
  {
    printf("\nCould not open file <%s>\n",filename);
    exit(0);
  }
  return 0;
}

// Allocate memory for matrix A[NROWS][MCOLS]
int alloc_2D(double *** A,double ** Adata, int NROWS, int MCOLS)
{
  int i;

  (*A)     = (double**) calloc(NROWS,sizeof(double*));
  (*Adata) = (double*) calloc(NROWS*MCOLS,sizeof(double)); 
  for (i = 0; i < NROWS; i++) (*A)[i] = (*Adata) + MCOLS*i;

  if ((*A) == NULL || (*Adata) == NULL)
  {
    printf("\nCould not allocate memory for 2D Array\n");
    exit(0);
  }

  return 0;
}

// Read a file containing matrix A[NROWS][MCOLS]
int read_file(FILE ** fp, double ***A, int NROWS, int MCOLS)
{
  int i,j;

  for (i = 0; i < (NROWS); i++)
    for (j = 0; j < (MCOLS); j++)
      fscanf((*fp),"%lf", (**A)+i*(MCOLS)+j);

  return 0;
}

// Print 2D Array A[NROWS][MCOLS]
int prn_2D(double ***A, int NROWS, int MCOLS)
{
  int i,j;

  for (i = 0; i < (NROWS); i++)
  {
    for (j = 0; j < (MCOLS); j++)
    {
      printf(" %10.5f", *((**A)+i*(MCOLS)+j));
      if (j == (MCOLS)-1)
        putchar('\n');
    }
  }

  return 0;
}

// Print 2D Array C[NROWS][MCOLS] to file
int fprn_2D(double ***C, int NROWS, int MCOLS, FILE ** fp, char filename[])
{
  int i,j;

  if (((*fp) = fopen(filename,"w")) == NULL){
    printf("\nError opening file <%s>.",filename);
    exit(0);
  }

  fprintf((*fp),"%d %d\n", NROWS, MCOLS);

  for (i = 0; i < NROWS; i++)
  {
    for (j = 0; j < MCOLS; j++)
    {
      fprintf((*fp)," %12.5e", *((**C)+i*(MCOLS)+j));
      if (j == (MCOLS-1))
        fputc('\n',(*fp));
    }
  }
  return 0;
}

// Set Boundary Conditions for T[i][j]
int setBC(double *** T, int NI, int NJ, double dx, double dy)
{
  int i,j;
  double x,y;

  for (i = 0; i < NI; i++)
  {
    x = (double) i*dx;
    for (j = 0; j < NJ; j++)
    {
      y = (double) j*dy;
      (*T)[0][j]    = 0.0;        // T(0,y)
      (*T)[NI-1][j] = 2.0*exp(y); // T(2,y)
      (*T)[i][0]    = x;          // T(x,0)
      (*T)[i][NJ-1] = exp(1.0) * x;    // T(x,1)
    }
  }
  return 0;
}

// Set Source Term S[i][j]
int setSourceTerm(double *** S, int NI, int NJ, double dx, double dy)
{
  int i,j;
  double x,y;

  for (i = 0; i < NI; i++)
  {
    x = (double) i*dx;
    for (j = 0; j < NJ; j++)
    {
      y = (double) j*dy;
      (*S)[i][j] = x*exp(y); // S(x,y)
    }
  }
  return 0;
}



// Compute inf-norm of ||Te-Tc||_inf
double inf_norm(double **T, double **S, int NI, int NJ)
{
  int i,j;
  double dt_max = 0.0;

  for (i = 1; i < NI-1; i++)
    for (j = 1; j < NJ-1; j++)
      dt_max = max(fabs(T[i][j]-S[i][j]),dt_max);

  return (dt_max);
}

/* ====== Print Results  ======= */
int prn_results( int NI, int NJ, double omega, double dt_max, int k_final, double delta_time)
{
  printf("==================================================\n");
  printf("  NI    NJ     w        E           k      time\n");
  printf("--------------------------------------------------\n");
  printf(" %3d", NI);
  printf(" %5d", NJ);
  printf(" %6.3f", omega);
  printf("  %12.5e", dt_max);
  printf(" %5d", k_final);
  printf(" %12lf\n\n",delta_time);

  return 0;
}

// print results to file
int prn_results2file( FILE ** fp, char filename[], int NI, int NJ,
           double dx, double dy, double omega, double dt_max, int k_final, double delta_time)
{
  if (((*fp) = fopen(filename,"w")) == NULL){
    printf("\nError opening file <%s>.",filename);
    exit(0);
  }
  fprintf((*fp),"NI    =  %3d\n", NI);
  fprintf((*fp),"NJ    =  %3d\n", NJ);
  fprintf((*fp),"dx    =  %.15e\n", dx);
  fprintf((*fp),"dy    =  %.15e\n", dy);
  fprintf((*fp),"error =  %.15e\n", dt_max);
  fprintf((*fp),"omega =  %6.3f\n", omega);
  fprintf((*fp),"k     =  %5d\n", k_final);
  fprintf((*fp),"time  = %12lf\n",delta_time);
  return 0;
}


// Write GNUPLOT file connectivity for triangle grid connectivity, nt = num triangles
int write_gnuplot(char buff[],char filename[],FILE *fp,int **tri_conn,int nt,double *x,double *y)
{
  int i;
  buff[0]='\0';
  strcat(buff,filename);
  strcat(buff,"_plot.dat");
  printf("\nGNUPLOT Filename = <%s>\n",buff);
  // Open file for write
  if ((fp = fopen(buff,"w")) == 0)
  {
    printf("\nError opening file <%s>.",buff);
    exit(0);
  }
  for (i=0; i < nt; i++)
  {
    int n0 = tri_conn[i][0];
    int n1 = tri_conn[i][1];
    int n2 = tri_conn[i][2];
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n0],y[n0]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n1],y[n1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n2],y[n2]);
    fprintf(fp,"%19.10e %19.10e 0.0\n\n",x[n0],y[n0]);
  }
  fclose(fp);
  fp = NULL;

  return 0;
}

// ################# Begin Function ###########################
// ################# Begin Function ###########################
int parse_args(char argc, char * argv[],int * maxiter, int * ask, int * ProblemType,
               double * M, double * alpha, char fname[], char ofile[])
{
  int i;
  int print = 0;
  int norm_flag = 0;
  int mach_flag = 0;
  int angle_flag = 0;
  int file_flag = 0;
  int ofile_flag = 0;
  int help_flag = 0;
  printf("number of args = %d\n", argc);

  // iterate over arguments
  for (i = 1; i < (argc); i++)
  {
    printf("argv[%d] = %s\n",i,argv[i]);
    if (strcmp("-h", argv[i]) == 0) {
      help_flag = 1;
      continue;
    }
    if (strcmp("-n", argv[i]) == 0){
      (*maxiter) = atoi(argv[++i]);
      norm_flag = 1;
      continue;
    }
    if (strcmp("-a", argv[i]) == 0){
      (*ask) = atoi(argv[++i]);
      norm_flag = 1;
      continue;
    }
    if (strcmp("-mach", argv[i]) == 0){
      (*M) = atof(argv[++i]);
      mach_flag = 1;
      continue;
    }
    if (strcmp("-alpha", argv[i]) == 0){
      (*alpha) = atof(argv[++i]);
      angle_flag = 1;
      continue;
    }
    if (strcmp("-f", argv[i]) == 0){
      strcat(fname,argv[++i]);
      file_flag = 1;
      continue;
    }
    if (strcmp("-o", argv[i]) == 0){
      strcat(ofile,argv[++i]);
      ofile_flag = 1;
      continue;
    }
    if (strcmp("-p", argv[i]) == 0) {
      (*ProblemType) = atoi(argv[++i]);
      print = 1; // Print T to filename
      continue;
    }
  }

#if 0
  // Get inputs that werent sent by command line
  if ((mach_flag < 1) || (angle_flag < 1))
  {
    printf("Enter Mach and alpha\n");
    scanf("%lf %lf", M, alpha);
  }
  if (norm_flag < 1)
  {
    printf("Enter xnorm and ynorm\n");
    scanf("%lf %lf", nx, ny);
  }
#endif
  if (help_flag || argc < 2)
  {
    printf("Quasi-1D Euler solver\n");
    printf("Options:\n");
    printf("n = num iters\n");
    printf("a = ask to continue\n");
    printf("f = input  grid file\n");
    printf("o = output grid file\n");
    //printf("P = Print\n");
    printf("h = help\n");
  }

  if (file_flag < 1)
  {
    printf("Enter 2D grid input filename\n");
    scanf("%s", fname);
  }

  if (print < 1)
  {
    while(*ProblemType > 2 || *ProblemType < 1)
    {
      printf("Enter problem type\n");
      printf("1: subsonic-supersonic\n");
      printf("2: subsonic\n");
      //printf("3: shock\n");
      scanf("%d", ProblemType);
    }
  }

  if (ofile_flag < 1)
    strcat(ofile,"csvsoln.csv");

  return 0;
}


//prn_tri(tri_conn,nt);
int prn_tri(int **tri_conn, int nt)
{
  int t;

  for (t=0; t < nt; t++)
  {
    printf("tri_conn[%3d][%3d %3d %3d]\n",t,tri_conn[t][0], tri_conn[t][1],tri_conn[t][2]);
  }
  return 0;
}


// ################# Begin Function ###########################
// ################# Begin Function ###########################
int prn_jac(double ** DFP, double ** DFM)
{

  printf("DFP\n");
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFP[0][0], DFP[0][1], DFP[0][2], DFP[0][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFP[1][0], DFP[1][1], DFP[1][2], DFP[1][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFP[2][0], DFP[2][1], DFP[2][2], DFP[2][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFP[3][0], DFP[3][1], DFP[3][2], DFP[3][3]);


  printf("DFM\n");
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFM[0][0], DFM[0][1], DFM[0][2], DFM[0][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFM[1][0], DFM[1][1], DFM[1][2], DFM[1][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFM[2][0], DFM[2][1], DFM[2][2], DFM[2][3]);
  printf(" % .7f  % .7f  % .7f  % .7f\n", DFM[3][0], DFM[3][1], DFM[3][2], DFM[3][3]);

  return 0;
}

// ################# Begin Function ###########################
// ################# Begin Function ###########################
int print_fluxes(double fp1, double fp2, double fp3, double fp4,
                 double fm1, double fm2, double fm3, double fm4)
{

  printf("F+\n");
  printf("----------------\n");
  printf("fplus(1)  = %18.15f\n", fp1);
  printf("fplus(2)  = %18.15f\n", fp2);
  printf("fplus(3)  = %18.15f\n", fp3);
  printf("fplus(4)  = %18.15f\n", fp4);

  printf("\nF-\n");
  printf("----------------\n");
  printf("fminus(1) = %18.15f\n", fm1);
  printf("fminus(2) = %18.15f\n", fm2);
  printf("fminus(3) = %18.15f\n", fm3);
  printf("fminus(4) = %18.15f\n", fm4);

  return 0;
}


