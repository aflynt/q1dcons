#include "flux.h"


// Calculate Q vector
// Q1,Q2,Q3,Q4 == Q_infinity conditions
int initQ(double *Q1,double *Q2,double *Q3,double *Q4, const double M, double *alpha)
{
  double rho, P, u, v, e;

  // degrees -> radians
  (*alpha) *= M_PI / 180.0;

  // Primitive Variables
  rho = 1.0;
  P = 1.0 / g;
  u = M*cos((*alpha));
  v = M*sin((*alpha));
  e = P / (g - 1.0) + 0.5*rho*(u*u + v*v);

  // Q vector
  (*Q1) = rho;
  (*Q2) = rho * u;
  (*Q3) = rho * v;
  (*Q4) = e;

  return 0;
}

int comp_resid(int nn, int nt, int nb, int **tri_conn, double * x, double *y,
               int *nbs, int ***bs, double **Q, double **R,
               double Q1, double Q2, double Q3, double Q4)
{
  //
  // Function updates R[node][4]
  // Q1,Q2,Q3,Q4 == Q_infinity conditions
  //
  int i,k,t,b;
  int    n0, n1, n2;
  double x0, x1, x2;
  double y0, y1, y2;
  double xc, yc, xm, ym;
  double fp1,fp2,fp3,fp4;
  double fm1,fm2,fm3,fm4;
  double nx = 0.0;
  double ny = 0.0;

  // loop over all triangles
  for (t = 0; t < nt; t++)
  {
    // loop over 3 nodes of triangle
    for (k = 0; k < 3; k++)
    {
      n0 = tri_conn[t][k];
      n1 = tri_conn[t][(k+1)%3];
      n2 = tri_conn[t][(k+2)%3];

      x0 = x[n0];
      x1 = x[n1];
      x2 = x[n2];

      y0 = y[n0];
      y1 = y[n1];
      y2 = y[n2];

      xc = (x0 + x1 + x2) / 3.0;
      yc = (y0 + y1 + y2) / 3.0;

      xm = 0.5*(x0 + x1);
      ym = 0.5*(y0 + y1);

      nx = yc - ym;
      ny = xm - xc;

      // get F+ and F-
      fluxp(Q[n0][0], Q[n0][1], Q[n0][2], Q[n0][3], nx, ny, &fp1, &fp2, &fp3, &fp4);
      fluxm(Q[n1][0], Q[n1][1], Q[n1][2], Q[n1][3], nx, ny, &fm1, &fm2, &fm3, &fm4);


      // Resid[Q_left]
      R[n0][0] += (fp1 + fm1);
      R[n0][1] += (fp2 + fm2);
      R[n0][2] += (fp3 + fm3);
      R[n0][3] += (fp4 + fm4);

      // Resid[Q_right]
      R[n1][0] -= (fp1 + fm1);
      R[n1][1] -= (fp2 + fm2);
      R[n1][2] -= (fp3 + fm3);
      R[n1][3] -= (fp4 + fm4);
    }
  }
  //printf("\nafter going over all triangles\n");
  //printf("R[node %3d][% 7.2f][% 7.2f][% 7.2f][% 7.2f]\n", 4, R[4][0], R[4][1], R[4][2], R[4][3]);

  // Loop over all boundaries
  for (b=0; b < nb; b++)
  {
    // Loop over boundary segments
    for (i=0; i < nbs[b]; i++)
    {
      //get node 0 1
      n0 = bs[b][i][0];
      n1 = bs[b][i][1];

      x0 = x[n0];
      x1 = x[n1];
      y0 = y[n0];
      y1 = y[n1];

      // get midpt
      xm = 0.5*(x0 + x1);
      ym = 0.5*(y0 + y1);

      // get normal from midpt to node 0
      nx = ym - y0;
      ny = x0 - xm;

      // fminus(Q_inf) for node 0 & 1
      fluxm(Q1,Q2,Q3,Q4,nx, ny, &fm1, &fm2, &fm3, &fm4);

      // fplus(QL) for node 0
      fluxp(0.75*Q[n0][0]+0.25*Q[n1][0],
            0.75*Q[n0][1]+0.25*Q[n1][1],
            0.75*Q[n0][2]+0.25*Q[n1][2],
            0.75*Q[n0][3]+0.25*Q[n1][3],
            nx, ny, &fp1, &fp2, &fp3, &fp4);

      // Update Resid for node 0
      R[n0][0] += (fp1 + fm1);
      R[n0][1] += (fp2 + fm2);
      R[n0][2] += (fp3 + fm3);
      R[n0][3] += (fp4 + fm4);


      // fplus(QL) for node 1
      fluxp(0.75*Q[n1][0]+0.25*Q[n0][0],
            0.75*Q[n1][1]+0.25*Q[n0][1],
            0.75*Q[n1][2]+0.25*Q[n0][2],
            0.75*Q[n1][3]+0.25*Q[n0][3],
            nx, ny, &fp1, &fp2, &fp3, &fp4);

      // Update Resid for node 1
      R[n1][0] += (fp1 + fm1);
      R[n1][1] += (fp2 + fm2);
      R[n1][2] += (fp3 + fm3);
      R[n1][3] += (fp4 + fm4);

    }
  }

  //printf("\nafter looping over boundaries\n");
  //printf("R[node %3d][% 7.2f][% 7.2f][% 7.2f][% 7.2f]\n", 4, R[4][0], R[4][1], R[4][2], R[4][3]);

  return 0;
}

