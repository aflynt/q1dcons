#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include "List.h"
//#define NDEBUG

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define dprint(x) printf(#x " = %3d\n",x);
#define ddprint(x) printf(#x " = % .7f\n",(double) x);
#define MIN_REQUIRED 2
#define CHECKPT {printf("Checkpoint: %s, line %d\n",__FILE__,__LINE__);\
fflush(stdout);}

// global gamma
#define g   1.4
#define gm1 0.4


int getn2n(int nn, int ntri, int **tri, List **nhash);
//int getCRS(int nn, List **nhash, int *IA, int *JA, int *IAU);
int getCRS(int nn, List **nhash, int **IA, int **JA, int **IAU);

int delta_t(int nn, int nt, int nb, int **tri_conn, double * x, double *y,
           int *nbs, int ***bs, double **Q, double *area, double *cdt);

int tri_area(int nn, int nt, double *x, double *y, int **tri_conn, double *area);

int comp_resid(int nn, int nt, int nb, int **tri_conn, double * x, double *y,
               int *nbs, int ***bs, double **Q, double **R,
               double Q1, double Q2, double Q3, double Q4);

int initQ(double *Q1,double *Q2,double *Q3,double *Q4, const double M, double *alpha);

int normal_flux(double   Q1, double   Q2, double   Q3, double   Q4, double nx, double ny,
                double *fp1, double *fp2, double *fp3, double *fp4);

int fluxp(double   Q1, double   Q2, double   Q3, double   Q4, double nx, double ny,
          double *fp1, double *fp2, double *fp3, double *fp4);
int fluxm(double   Q1, double   Q2, double   Q3, double   Q4, double nx, double ny,
          double *fm1, double *fm2, double *fm3, double *fm4);

int ffp(double   Q1, double   Q2, double   Q3, double   Q4, double nx, double ny,
          double *fp1, double *fp2, double *fp3, double *fp4);
int ffm(double Q1, double Q2, double Q3, double Q4, double nxh, double nyh,
          double *fm1, double *fm2, double *fm3, double *fm4);


int jacp(double Q1, double Q2, double Q3,  double Q4,
         double nx, double ny, double fp1, double ** DFP);
int jacm(double Q1, double Q2, double Q3,  double Q4,
         double nx, double ny, double fm1, double ** DFM);

int print_fluxes(double fp1, double fp2, double fp3, double fp4,
                 double fm1, double fm2, double fm3, double fm4);

int prn_jac(double ** DFP, double ** DFM);

int prn_tri(int **tri_conn, int nt);

int parse_args(char argc, char * argv[],int * maxiter, int * ask,
               double * M, double * alpha, char fname[], char ofile[]);


int help(void);
int open_file(char * filename, FILE ** fp);
int alloc_2D(double *** A,double ** Adata, int NROWS, int MCOLS);
int read_file(FILE ** fp, double ***A, int NROWS, int MCOLS);
int prn_2D(double ***A, int NROWS, int MCOLS);
int fprn_2D(double ***C, int NROWS, int MCOLS, FILE ** fp, char filename[]);
int setBC(double *** T, int NI, int NJ, double dx, double dy);
int setSourceTerm(double *** S, int NI, int NJ, double dx, double dy);
double inf_norm(double **T, double **S, int NI, int NJ);
int prn_results( int NI, int NJ, double omega, double dt_max, int k_final, double delta_time);
int prn_results2file( FILE ** fp, char filename[], int NI, int NJ,
           double dx, double dy, double omega, double dt_max, int k_final, double delta_time);
int write_gnuplot(char buff[],char filename[],FILE *fp,int **tri_conn,int nt,double *x,double *y);
