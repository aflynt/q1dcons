#include "flux.h"



/* === Triangle Area === */
// area[nn] == node area
int tri_area(int nn, int nt, double *x, double *y, int **tri_conn, double *area)
{
  double a, ux, uy, vx, vy;
  int n, t, n0, n1, n2;

  double amn =  1.0e10;
  double amx = -1.0e10;
  double avg =  0.0;

  for (n=0; n < nn; n++) // zero out node_areas
    area[n] = 0.0;

  for (t=0; t < nt; t++) // loop over all triangles
  {
    n0 = tri_conn[t][0];
    n1 = tri_conn[t][1];
    n2 = tri_conn[t][2];
    ux = x[n1]-x[n0];
    uy = y[n1]-y[n0];
    vx = x[n2]-x[n0];
    vy = y[n2]-y[n0];

    // get area of triangle
    a = 0.5*(ux*vy-uy*vx);

    // add triangles' area to nodes
    area[n0] += a / 3.0;
    area[n1] += a / 3.0;
    area[n2] += a / 3.0;

    //printf("tri: %2d Area = %.6f\n", t, a);

    amn = min(a,amn);
    amx = max(a,amx);
    avg += a;
  }
#if 0
  printf("\nTriangle minimum area = %17.10e",amn);
  printf("\nTriangle average area = %17.10e",avg/max(1,nt));
  printf("\nTriangle maximum area = %17.10e\n",amx);
#endif

  return 0;
}




//============================== delta_t =====================================
//
//  Calculate a time step for each node
//  Q is the conservative variable
//
//============================================================================
int delta_t(int nn, int nt, int nb, int **tri_conn, double * x, double *y,
           int *nbs, int ***bs, double **Q, double *area, double *cdt)
{
  int i,k,n,t,b;
  int    n0, n1, n2;
  double x0, x1, x2;
  double y0, y1, y2;
  double xc, yc, xm, ym;
  double nx = 0.0;
  double ny = 0.0;
  double Vn;

  double nxh,nyh,l;
  double u,v,c;
  double rho0,e0,c0,u0,v0,P0;
  double rho1,e1,c1,u1,v1,P1;
  double ubar;


  //Zero out cdt
  for(n=0; n<nn; n++)
    cdt[n] = 0.0;


  // loop over all triangles
  for (t = 0; t < nt; t++)
  {
    // loop over 3 nodes of triangle
    for (k = 0; k < 3; k++)
    {
      // get nodes of triangle
      n0 = tri_conn[t][k];
      n1 = tri_conn[t][(k+1)%3];
      n2 = tri_conn[t][(k+2)%3];

      // get node coordinates
      x0 = x[n0];
      x1 = x[n1];
      x2 = x[n2];

      y0 = y[n0];
      y1 = y[n1];
      y2 = y[n2];

      // get center coords
      xc = (x0 + x1 + x2) / 3.0;
      yc = (y0 + y1 + y2) / 3.0;

      // get edge midpt
      xm = 0.5*(x0 + x1);
      ym = 0.5*(y0 + y1);

      // calculate normals
      nx = yc - ym;
      ny = xm - xc;
      l = sqrt(nx*nx + ny*ny);
      nxh = nx/l;
      nyh = ny/l;

      // Flow variables
      rho0 = Q[n0][0];
      u0   = Q[n0][1] / rho0;
      v0   = Q[n0][2] / rho0;
      e0   = Q[n0][3];
      P0 = gm1*(e0 - 0.5*rho0*(u0*u0 + v0*v0));
      c0 = sqrt(g*P0/rho0);

      rho1 = Q[n1][0];
      u1   = Q[n1][1] / rho1;
      v1   = Q[n1][2] / rho1;
      e1   = Q[n1][3];
      P1 = gm1*(e1 - 0.5*rho1*(u1*u1 + v1*v1));
      c1 = sqrt(g*P1/rho1);

      u = 0.5*(u0 + u1);
      v = 0.5*(v0 + v1);
      c = 0.5*(c0 + c1);

      ubar = u*nxh + v*nyh;

      Vn = (fabs(ubar) + c) * l;

      // add to time step for nodes
      cdt[n0] += Vn;
      cdt[n1] += Vn;

    }
  }

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

      // Flow quantities
      rho0 = Q[n0][0];
      u0   = Q[n0][1] / rho0;
      v0   = Q[n0][2] / rho0;
      e0   = Q[n0][3];
      P0 = gm1*(e0 - 0.5*rho0*(u0*u0 + v0*v0));
      c0 = sqrt(g*P0/rho0);

      rho1 = Q[n1][0];
      u1   = Q[n1][1] / rho1;
      v1   = Q[n1][2] / rho1;
      e1   = Q[n1][3];
      P1 = gm1*(e1 - 0.5*rho1*(u1*u1 + v1*v1));
      c1 = sqrt(g*P1/rho1);

      // Left side of face
      nx =   ym - y0;
      ny = -(xm - x0);
      l = sqrt(nx*nx + ny*ny);
      u = 5./6.*u0 + 1./6.*u1;
      v = 5./6.*v0 + 1./6.*v1;
      c = 5./6.*c0 + 1./6.*c1;

      Vn = fabs(u*nx + v*ny) + c*l;

      cdt[n0] += Vn;

      // Right side of face
      nx =   y1 - ym;
      ny = -(x1 - xm);
      l = sqrt(nx*nx + ny*ny);
      u = 5./6.*u1 + 1./6.*u0;
      v = 5./6.*v1 + 1./6.*v0;
      c = 5./6.*c1 + 1./6.*c0;

      Vn = fabs(u*nx + v*ny) + c*l;

      cdt[n1] += Vn;
    }
  }

  // cdt = area / cdt
  for(n=0; n<nn; n++)
    cdt[n] = area[n] / cdt[n];

  return 0;
}
