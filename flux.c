#include "flux.h"


int fluxp(double Q1, double Q2, double Q3, double Q4, double nx, double ny,
          double *fp1, double *fp2, double *fp3, double *fp4)
{
  double nxh,nyh,l;
  double rho;
  double e,c,u,v,P;
  double ubar;
  double M;

  //1 calculate primitive variables
  rho = Q1;
  e = Q4;
  u = Q2/Q1;
  v = Q3/Q1;

  //2 calculate normals
  l = sqrt(nx*nx + ny*ny);
  nxh = nx/l;
  nyh = ny/l;

  //3 Pressure
  P = (g-1.0)*(e - 0.5*rho*(u*u + v*v));
  c = sqrt(g*P/rho);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  // M >= +1
  if ( M >= 1.0 )
  {
    //fplus = f
    (*fp1) = rho * ubar * l;
    (*fp2) = (rho*u*ubar + nxh*P) * l;
    (*fp3) = (rho*v*ubar + nyh*P) * l;
    (*fp4) = (e + P)*ubar * l;
  }

  // M <= -1
  else if ( M <= -1.0)
  {
    //fplus = 0
    (*fp1) = 0.0;
    (*fp2) = 0.0;
    (*fp3) = 0.0;
    (*fp4) = 0.0;
  }

  // |Mach| < 1
  else
  {
    // fluxp
    (*fp1) = 0.25*rho*c*(ubar/c+1.0)*(ubar/c+1.0)*l;
    (*fp2) = (*fp1)*(nxh/g*(-ubar + 2.0*c)+u);
    (*fp3) = (*fp1)*(nyh/g*(-ubar + 2.0*c)+v);
    (*fp4) = (*fp1)*((-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c)/(g*g-1.0) + 0.5*(u*u + v*v));
  }

  return 0;
}

// flux-
// ################# Begin Function ###########################
// ################# Begin Function ###########################
int fluxm(double Q1, double Q2, double Q3, double Q4, double nx, double ny,
          double *fm1, double *fm2, double *fm3, double *fm4)
{
  double nxh,nyh,l;
  double rho;
  double e,c,u,v,P;
  double ubar;
  double M;

  //1 calculate primitive variables
  rho = Q1;
  e = Q4;
  u = Q2/Q1;
  v = Q3/Q1;

  //2 calculate normals
  l = sqrt(nx*nx + ny*ny);
  nxh = nx/l;
  nyh = ny/l;

  //3 Pressure
  P = (g-1.0)*(e - 0.5*rho*(u*u + v*v));
  c = sqrt(g*P/rho);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  // M >= +1
  if ( M >= 1.0 )
  {
    //fminus = 0
    (*fm1) = 0.0;
    (*fm2) = 0.0;
    (*fm3) = 0.0;
    (*fm4) = 0.0;
  }

  // M <= -1
  else if ( M <= -1.0 )
  {
    //fminus = f
    (*fm1) = rho * ubar * l;
    (*fm2) = (rho*u*ubar + nxh*P) * l;
    (*fm3) = (rho*v*ubar + nyh*P) * l;
    (*fm4) = (e + P)*ubar * l;
  }

  // |Mach| < 1
  else
  {
    // fluxm
    (*fm1) = -0.25*rho*c*(ubar/c - 1.0)*(ubar/c - 1.0)*l;
    (*fm2) = (*fm1)*(nxh/g*(-ubar - 2.0*c)+u);
    (*fm3) = (*fm1)*(nyh/g*(-ubar - 2.0*c)+v);
    (*fm4) = (*fm1)*((-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c)/(g*g-1.0) + 0.5*(u*u + v*v));
  }

  return 0;
}

// jacobian function
// ################# Begin Function ###########################
// ################# Begin Function ###########################
int jacp(double Q1, double Q2, double Q3, double Q4,
         double nx, double ny, double fp1, double ** DFP)
{
  double nxh,nyh,l;
  double rho;
  double e,c,u,v,P;
  double ubar;
  double M;

  //1 calculate primitive variables
  rho = Q1;
  e = Q4;
  u = Q2/Q1;
  v = Q3/Q1;

  //2 calculate normals
  l = sqrt(nx*nx + ny*ny);
  nxh = nx/l;
  nyh = ny/l;

  //3 Pressure
  P = (g-1.0)*(e - 0.5*rho*(u*u + v*v));
  c = sqrt(g*P/rho);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  // Calc DFP[0][0]
  double A,B;
  double dA, dB;
  double u_Q1 = -Q2/(Q1*Q1);
  double v_Q1 = -Q3/(Q1*Q1);
  double c_Q1 = 0.25*g*(g-1.0)*(u*u+v*v)/(rho*c) - 0.5*c/rho;
  double ubar_Q1 = -(nxh*Q2/(Q1*Q1) + nyh*Q3/(Q1*Q1));
  double P_Q1 = 0.5*(g-1.0)*(u*u+v*v);

  double u_Q2 = 1.0 / Q1;
  double v_Q2 = 0.0;
  double c_Q2 = -g*(g-1.0)*Q2 / (2.0*c*Q1*Q1);
  double ubar_Q2 = nxh / Q1;
  double P_Q2 = -(g-1.0)*u;

  double u_Q3 = 0.0;
  double v_Q3 = 1.0 / Q1;
  double c_Q3 = -g*(g-1.0)*Q3 / (2.0*c*Q1*Q1);
  double ubar_Q3 = nyh / Q1;
  double P_Q3 = -(g-1.0)*v;

  double u_Q4 = 0.0;
  double v_Q4 = 0.0;
  double c_Q4 = g*(g-1.0) / (2.0*c*Q1*Q1);
  double ubar_Q4 = 0.0;


// DFP 1*
// =============== //

  // DFP 11
  A = 0.25*rho*c;
  B = (M+1.0)*(M+1.0);
  dA = 0.25*(rho*c_Q1 + c);
  dB = 2.0*(M+1.0)*(c*ubar_Q1-ubar*c_Q1)/(c*c);
  DFP[0][0] = (A*dB + B*dA)*l;

  // DFP 12
  A = 0.25*rho*c;
  B = (M+1.0)*(M+1.0);
  dA = 0.25*rho*c_Q2;
  dB = 2.0*(M+1.0)*(c*ubar_Q2 - ubar*c_Q2)/(c*c);
  DFP[0][1] = (A*dB + B*dA)*l;

  // DFP 13
  A = 0.25*rho*c;
  B = (M+1.0)*(M+1.0);
  dA = 0.25*rho*c_Q3;
  dB = 2.0*(M+1.0)*(c*ubar_Q3 - ubar*c_Q3)/(c*c);
  DFP[0][2] = (A*dB + B*dA)*l;

  // DFP 14
  A = 0.25*rho*c;
  B = (M+1.0)*(M+1.0);
  dA = 0.25*rho*c_Q4;
  dB = 2.0*(M+1.0)*(-ubar*c_Q4)/(c*c);
  DFP[0][3] = (A*dB + B*dA)*l;

// DFP 2*
// =============== //

  // DFP 21
  B = (nxh/g*(-ubar + 2.0*c)+u);
  dB = nxh/g*(-ubar_Q1 + 2.0*c_Q1) + u_Q1;
  DFP[1][0] = fp1*dB + DFP[0][0]*B;

  // DFP 22
  B = (nxh/g*(-ubar + 2.0*c)+u);
  dB = nxh/g*(-ubar_Q2 + 2.0*c_Q2) + u_Q2;
  DFP[1][1] = fp1*dB + DFP[0][1]*B;

  // DFP 23
  B = (nxh/g*(-ubar + 2.0*c)+u);
  dB = nxh/g*(-ubar_Q3 + 2.0*c_Q3) + u_Q3;
  DFP[1][2] = fp1*dB + DFP[0][2]*B;

  // DFP 24
  B = (nxh/g*(-ubar + 2.0*c)+u);
  dB = nxh/g*2.0*c_Q4;
  DFP[1][3] = fp1*dB + DFP[0][3]*B;


// DFP 3*
// =============== //

  // DFP 31
  B = (nyh/g*(-ubar + 2.0*c)+v);
  dB = nyh/g*(-ubar_Q1 + 2.0*c_Q1) + v_Q1;
  DFP[2][0] = fp1*dB + DFP[0][0]*B;

  // DFP 32
  B = (nyh/g*(-ubar + 2.0*c)+v);
  dB = nyh/g*(-ubar_Q2 + 2.0*c_Q2) + v_Q2;
  DFP[2][1] = fp1*dB + DFP[0][1]*B;

  // DFP 33
  B = (nyh/g*(-ubar + 2.0*c)+v);
  dB = nyh/g*(-ubar_Q3 + 2.0*c_Q3) + v_Q3;
  DFP[2][2] = fp1*dB + DFP[0][2]*B;

  // DFP 34
  B = (nyh/g*(-ubar + 2.0*c)+v);
  dB = nyh/g*2.0*c_Q4;
  DFP[2][3] = fp1*dB + DFP[0][3]*B;


// DFP 4*
// =============== //

  // DFP 41
  B = (-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q1 + 2.0*(g-1.0)*(ubar*c_Q1 + c*ubar_Q1) + 4.0*c*c_Q1) / (g*g - 1.0)
          + u*u_Q1 + v*v_Q1;
  DFP[3][0] = fp1*dB + DFP[0][0]*B;

  // DFP 42
  B = (-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q2 + 2.0*(g-1.0)*(ubar*c_Q2 + c*ubar_Q2) + 4.0*c*c_Q2) / (g*g - 1.0)
          + u*u_Q2 + v*v_Q2;
  DFP[3][1] = fp1*dB + DFP[0][1]*B;

  // DFP 43
  B = (-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q3 + 2.0*(g-1.0)*(ubar*c_Q3 + c*ubar_Q3) + 4.0*c*c_Q3) / (g*g - 1.0)
          + u*u_Q3 + v*v_Q3;
  DFP[3][2] = fp1*dB + DFP[0][2]*B;

  // DFP 44
  B = (-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (2.0*(g-1.0)*ubar*c_Q4 + 4.0*c*c_Q4) / (g*g - 1.0);
  DFP[3][3] = fp1*dB + DFP[0][3]*B;

  // Normal flux calculations
  if (M >= 1.0)
  {
    DFP[0][0] = (rho*ubar_Q1 + ubar)*l;
    DFP[0][1] = (rho*ubar_Q2)*l;
    DFP[0][2] = (rho*ubar_Q3)*l;
    DFP[0][3] = (rho*ubar_Q4)*l;

    DFP[1][0] = (ubar*(u + u_Q1*rho) + ubar_Q1*rho*u + nxh*P_Q1)*l;
    DFP[1][1] = (ubar*(u_Q2*rho) + ubar_Q2*rho*u + nxh*P_Q2)*l;
    DFP[1][2] = (ubar*(u_Q3*rho) + ubar_Q3*rho*u + nxh*P_Q3)*l;
    DFP[1][3] = nxh*(g-1.0)*l;

    DFP[2][0] = ( rho*(v*ubar_Q1 + ubar*v_Q1) + v*ubar + nyh*P_Q1 )*l;
    DFP[2][1] = ( rho*(v*ubar_Q2) + nyh*P_Q2 )*l;
    DFP[2][2] = ( rho*(v*ubar_Q3 + ubar*v_Q3) + nyh*P_Q3 )*l;
    DFP[2][3] = nyh*(g-1.0)*l;

    DFP[3][0] = ( (e+P)*ubar_Q1 + ubar*P_Q1 )*l;
    DFP[3][1] = ( (e+P)*ubar_Q2 + ubar*P_Q2 )*l;
    DFP[3][2] = ( (e+P)*ubar_Q3 + ubar*P_Q3 )*l;
    DFP[3][3] = ubar*g*l;
  }


  // M less than -1.0
  else if (M <= -1.0)
  {
    DFP[0][0] = 0.0;
    DFP[0][1] = 0.0;
    DFP[0][2] = 0.0;
    DFP[0][3] = 0.0;

    DFP[1][0] = 0.0;
    DFP[1][1] = 0.0;
    DFP[1][2] = 0.0;
    DFP[1][3] = 0.0;

    DFP[2][0] = 0.0;
    DFP[2][1] = 0.0;
    DFP[2][2] = 0.0;
    DFP[2][3] = 0.0;

    DFP[3][0] = 0.0;
    DFP[3][1] = 0.0;
    DFP[3][2] = 0.0;
    DFP[3][3] = 0.0;
  }

  return 0;
}

// jacobian function
// ################# Begin Function ###########################
// ################# Begin Function ###########################
int jacm(double Q1, double Q2, double Q3, double Q4,
             double nx, double ny, double fm1, double ** DFM)
{
  double nxh,nyh,l;
  double rho;
  double e,c,u,v,P;
  double ubar;
  double M;

  //1 calculate primitive variables
  rho = Q1;
  e = Q4;
  u = Q2/Q1;
  v = Q3/Q1;

  //2 calculate normals
  l = sqrt(nx*nx + ny*ny);
  nxh = nx/l;
  nyh = ny/l;

  //3 Pressure
  P = (g-1.0)*(e - 0.5*rho*(u*u + v*v));
  c = sqrt(g*P/rho);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  double A,B;
  double dA, dB;
  double u_Q1 = -Q2/(Q1*Q1);
  double v_Q1 = -Q3/(Q1*Q1);
  double c_Q1 = 0.25*g*(g-1.0)*(u*u+v*v)/(rho*c) - 0.5*c/rho;
  double ubar_Q1 = -(nxh*Q2/(Q1*Q1) + nyh*Q3/(Q1*Q1));
  double P_Q1 = 0.5*(g-1.0)*(u*u+v*v);

  double u_Q2 = 1.0 / Q1;
  double v_Q2 = 0.0;
  double c_Q2 = -g*(g-1.0)*Q2 / (2.0*c*Q1*Q1);
  double ubar_Q2 = nxh / Q1;
  double P_Q2 = -(g-1.0)*u;

  double u_Q3 = 0.0;
  double v_Q3 = 1.0 / Q1;
  double c_Q3 = -g*(g-1.0)*Q3 / (2.0*c*Q1*Q1);
  double ubar_Q3 = nyh / Q1;
  double P_Q3 = -(g-1.0)*v;

  double u_Q4 = 0.0;
  double v_Q4 = 0.0;
  double c_Q4 = g*(g-1.0) / (2.0*c*Q1*Q1);
  double ubar_Q4 = 0.0;


// DF 1*
//############################################

  // DFM 11
  A = -0.25*rho*c;
  B = (M-1.0)*(M-1.0);
  dA = -0.25*(rho*c_Q1 + c);
  dB = 2.0*(M-1.0)*(c*ubar_Q1 - ubar*c_Q1)/(c*c);
  DFM[0][0] = (A*dB + B*dA)*l;

  // DFM 12
  A = -0.25*rho*c;
  B = (M-1.0)*(M-1.0);
  dA = -0.25*rho*c_Q2;
  dB = 2.0*(M-1.0)*(c*ubar_Q2 - ubar*c_Q2)/(c*c);
  DFM[0][1] = (A*dB + B*dA)*l;

  // DFM 13
  A = -0.25*rho*c;
  B = (M-1.0)*(M-1.0);
  dA = -0.25*rho*c_Q3;
  dB = 2.0*(M-1.0)*(c*ubar_Q3 - ubar*c_Q3)/(c*c);
  DFM[0][2] = (A*dB + B*dA)*l;

  // DFM 14
  A = -0.25*rho*c;
  B = (M-1.0)*(M-1.0);
  dA = -0.25*rho*c_Q4;
  dB = 2.0*(M-1.0)*(-ubar*c_Q4)/(c*c);
  DFM[0][3] = (A*dB + B*dA)*l;


// DF 2*
//############################################

  // DFM 21
  B = (nxh/g*(-ubar - 2.0*c)+u);
  dB = nxh/g*(-ubar_Q1 - 2.0*c_Q1) + u_Q1;
  DFM[1][0] = fm1*dB + DFM[0][0]*B;

  // DFM 22
  B = (nxh/g*(-ubar - 2.0*c)+u);
  dB = nxh/g*(-ubar_Q2 - 2.0*c_Q2) + u_Q2;
  DFM[1][1] = fm1*dB + DFM[0][1]*B;

  // DFM 23
  B = (nxh/g*(-ubar - 2.0*c)+u);
  dB = nxh/g*(-ubar_Q3 - 2.0*c_Q3) + u_Q3;
  DFM[1][2] = fm1*dB + DFM[0][2]*B;

  // DFM 24
  B = (nxh/g*(-ubar - 2.0*c)+u);
  dB = nxh/g*(-2.0*c_Q4);
  DFM[1][3] = fm1*dB + DFM[0][3]*B;


// DF 3*
//############################################

  // DFM 31
  B = (nyh/g*(-ubar - 2.0*c)+v);
  dB = nyh/g*(-ubar_Q1 - 2.0*c_Q1) + v_Q1;
  DFM[2][0] = fm1*dB + DFM[0][0]*B;

  // DFM 32
  B = (nyh/g*(-ubar - 2.0*c)+v);
  dB = nyh/g*(-ubar_Q2 - 2.0*c_Q2) + v_Q2;
  DFM[2][1] = fm1*dB + DFM[0][1]*B;

  // DFM 33
  B = (nyh/g*(-ubar - 2.0*c)+v);
  dB = nyh/g*(-ubar_Q3 - 2.0*c_Q3) + v_Q3;
  DFM[2][2] = fm1*dB + DFM[0][2]*B;

  // DFM 34
  B = (nyh/g*(-ubar - 2.0*c)+v);
  dB = nyh/g*(-2.0*c_Q4);
  DFM[2][3] = fm1*dB + DFM[0][3]*B;


// DF *4
//############################################

  // DFM 41
  B = (-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q1 - 2.0*(g-1.0)*(ubar*c_Q1 + c*ubar_Q1) + 4.0*c*c_Q1) / (g*g - 1.0)
          + u*u_Q1 + v*v_Q1;
  DFM[3][0] = fm1*dB + DFM[0][0]*B;

  // DFM 42
  B = (-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q2 - 2.0*(g-1.0)*(ubar*c_Q2 + c*ubar_Q2) + 4.0*c*c_Q2) / (g*g - 1.0)
          + u*u_Q2 + v*v_Q2;
  DFM[3][1] = fm1*dB + DFM[0][1]*B;

  // DFM 43
  B = (-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*ubar_Q3 - 2.0*(g-1.0)*(ubar*c_Q3 + c*ubar_Q3) + 4.0*c*c_Q3) / (g*g - 1.0)
          + u*u_Q3 + v*v_Q3;
  DFM[3][2] = fm1*dB + DFM[0][2]*B;

  // DFM 44
  B = (-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c) / (g*g - 1.0) + 0.5*(u*u + v*v);
  dB = (-2.0*(g-1.0)*ubar*c_Q4 + 4.0*c*c_Q4) / (g*g - 1.0);
  DFM[3][3] = fm1*dB + DFM[0][3]*B;


  // Normal flux calculations
  if (M >= 1.0)
  {
    DFM[0][0] = 0.0;
    DFM[0][1] = 0.0;
    DFM[0][2] = 0.0;
    DFM[0][3] = 0.0;

    DFM[1][0] = 0.0;
    DFM[1][1] = 0.0;
    DFM[1][2] = 0.0;
    DFM[1][3] = 0.0;

    DFM[2][0] = 0.0;
    DFM[2][1] = 0.0;
    DFM[2][2] = 0.0;
    DFM[2][3] = 0.0;

    DFM[3][0] = 0.0;
    DFM[3][1] = 0.0;
    DFM[3][2] = 0.0;
    DFM[3][3] = 0.0;
  }


  // M less than -1.0
  else if (M <= -1.0)
  {
    DFM[0][0] = (rho*ubar_Q1 + ubar)*l;
    DFM[0][1] = (rho*ubar_Q2)*l;
    DFM[0][2] = (rho*ubar_Q3)*l;
    DFM[0][3] = (rho*ubar_Q4)*l;

    DFM[1][0] = (ubar*(u + u_Q1*rho) + ubar_Q1*rho*u + nxh*P_Q1)*l;
    DFM[1][1] = (ubar*(u_Q2*rho) + ubar_Q2*rho*u + nxh*P_Q2)*l;
    DFM[1][2] = (ubar*(u_Q3*rho) + ubar_Q3*rho*u + nxh*P_Q3)*l;
    DFM[1][3] = nxh*(g-1.0)*l;

    DFM[2][0] = ( rho*(v*ubar_Q1 + ubar*v_Q1) + v*ubar + nyh*P_Q1 )*l;
    DFM[2][1] = ( rho*(v*ubar_Q2) + nyh*P_Q2 )*l;
    DFM[2][2] = ( rho*(v*ubar_Q3 + ubar*v_Q3) + nyh*P_Q3 )*l;
    DFM[2][3] = nyh*(g-1.0)*l;

    DFM[3][0] = ( (e+P)*ubar_Q1 + ubar*P_Q1 )*l;
    DFM[3][1] = ( (e+P)*ubar_Q2 + ubar*P_Q2 )*l;
    DFM[3][2] = ( (e+P)*ubar_Q3 + ubar*P_Q3 )*l;
    DFM[3][3] = ubar*g*l;
  }

  return 0;
}


// Alternative forms of functions. Investigated for speed
// flux+
// ################# Begin Function ###########################
// ################# Begin Function ###########################
int normal_flux(double Q1, double Q2, double Q3, double Q4, double nx, double ny,
                double *fp1, double *fp2, double *fp3, double *fp4)
{
  double nxh,nyh,l;
  double rho;
  double e,c,u,v,P;
  double ubar;
  double M;

  //1 calculate primitive variables
  rho = Q1;
  e = Q4;
  u = Q2/Q1;
  v = Q3/Q1;

  //2 calculate normals
  l = sqrt(nx*nx + ny*ny);
  nxh = nx/l;
  nyh = ny/l;

  //3 Pressure
  P = (g-1.0)*(e - 0.5*rho*(u*u + v*v));
  c = sqrt(g*P/rho);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  //fplus = f
  (*fp1) = rho * ubar * l;
  (*fp2) = (rho*u*ubar + nxh*P) * l;
  (*fp3) = (rho*v*ubar + nyh*P) * l;
  (*fp4) = (e + P)*ubar * l;

  return 0;
}

int ffp(double Q1, double Q2, double Q3, double Q4, double nxh, double nyh,
          double *fp1, double *fp2, double *fp3, double *fp4)
{
  double c,P;
  double ubar;
  double M;
  double u = Q2/Q1;
  double v = Q3/Q1;

  double l = sqrt(nxh*nxh + nyh*nyh);
  nxh /= l;
  nyh /= l;

  //3 Pressure
  P = (g-1.0)*(   Q4 - 0.5*Q1*(u*u + v*v)  );
  c = sqrt(g*P/Q1);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  // M >= +1
  if ( M >= 1.0 )
  {
    //fplus = f
    (*fp1) = Q1 * ubar * l;
    (*fp2) = (Q1*u*ubar + nxh*P) * l;
    (*fp3) = (Q1*v*ubar + nyh*P) * l;
    (*fp4) = (Q4 + P)*ubar * l;
  }

  // M <= -1
  else if ( M <= -1.0)
  {
    //fplus = 0
    (*fp1) = 0.0;
    (*fp2) = 0.0;
    (*fp3) = 0.0;
    (*fp4) = 0.0;
  }

  // |Mach| < 1
  else
  {
    // fluxp
    (*fp1) = 0.25*Q1*c*(ubar/c+1.0)*(ubar/c+1.0) * l;
    (*fp2) = (*fp1)*(nxh/g*(-ubar + 2.0*c)+u);
    (*fp3) = (*fp1)*(nyh/g*(-ubar + 2.0*c)+v);
    (*fp4) = (*fp1)*((-(g-1.0)*ubar*ubar + 2.0*(g-1.0)*ubar*c + 2.0*c*c)/(g*g-1.0) + 0.5*(u*u + v*v));
  }

  return 0;
}

int ffm(double Q1, double Q2, double Q3, double Q4, double nxh, double nyh,
          double *fm1, double *fm2, double *fm3, double *fm4)
{
  double c,P;
  double ubar;
  double M;
  double u = Q2/Q1;
  double v = Q3/Q1;

  double l = sqrt(nxh*nxh + nyh*nyh);
  nxh /= l;
  nyh /= l;

  //3 Pressure
  P = (g-1.0)*(Q4 - 0.5*Q1*(u*u + v*v));
  c = sqrt(g*P/Q1);

  //4 U_bar
  ubar = nxh*u + nyh*v;

  M = ubar / c;
  //ddprint(M);

  // M >= +1
  if ( M >= 1.0 )
  {
    //fminus = 0
    (*fm1) = 0.0;
    (*fm2) = 0.0;
    (*fm3) = 0.0;
    (*fm4) = 0.0;
  }

  // M <= -1
  else if ( M <= -1.0 )
  {
    //fminus = f
    (*fm1) = Q1 * ubar * l;
    (*fm2) = (Q1*u*ubar + nxh*P) * l;
    (*fm3) = (Q1*v*ubar + nyh*P) * l;
    (*fm4) = (Q4 + P)*ubar * l;
  }

  // |Mach| < 1
  else
  {
    // fluxm
    (*fm1) = -0.25*Q1*c*(ubar/c - 1.0)*(ubar/c - 1.0) * l;
    (*fm2) = (*fm1)*(nxh/g*(-ubar - 2.0*c)+u);
    (*fm3) = (*fm1)*(nyh/g*(-ubar - 2.0*c)+v);
    (*fm4) = (*fm1)*((-(g-1.0)*ubar*ubar - 2.0*(g-1.0)*ubar*c + 2.0*c*c)/(g*g-1.0) + 0.5*(u*u + v*v));
  }

  return 0;
}
