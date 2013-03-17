#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
   Framework code for reading 2-D generic unstructured mesh and performing
   metric calculations.
*/

double half(double p1, double p2);
double sq(double val);

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define PI 3.141592653589793

int main(int argcs, char* pArgs[])
{
  FILE *fp;
  int bdim = 80;
  char buff[bdim];
  int i, n, nn, nq, nt, q, t, nblocks;
  int **tri_conn, **quad_conn;
  double *x, *y;

  // Check for correct number of arguments (input and output file names)
  if (argcs != 2)
  {
    printf("\nUsage:");
    printf("\nmetrics input_2D_grid_file\n\n");
    exit(0);
  }

  /* Begin reading 2D grid file */

  /* Open 2D grid file for reading */
  if ((fp=fopen(pArgs[1],"r")) == NULL)
  {
    printf("\nCould not open file <%s>\n",pArgs[1]);
    exit(0);
  }

  /* Read number of grid points */
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of grid points
  sscanf(buff,"%d",&nn);
  printf("\nNumber of grid points = %d",nn);
  if ((x = (double*)malloc(nn*sizeof(double))) == NULL)
  {
    printf("\nCould not allocate memory for X");
    exit(0);
  }
  if ((y = (double*)malloc(nn*sizeof(double))) == NULL)
  {
    printf("\nCould not allocate memory for Y");
    exit(0);
  }
  /* Read 2D grid nodes */
  for (n=0; n < nn; n++)
  {
    fgets(buff,bdim,fp);
    sscanf(buff,"%lg %lg",&x[n],&y[n]);
  }
  
  /* Read number of blocks */
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of blocks
  sscanf(buff,"%d",&nblocks);
  if (nblocks != 1)
  {
    printf("\nNumber of blocks should be 1");
    exit(0);
  }

  /* Read number of triangular elements */
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of tri-elements
  sscanf(buff,"%d",&nt);
  printf("\nNumber of triangles = %d",nt);
  if ((tri_conn = (int**)malloc(nt*sizeof(int*))) == NULL)
  {
    printf("\nCould not allocate memory for triangle connectivity");
    exit(0);
  }
  for (t=0; t < nt; t++)
  {
    if ((tri_conn[t] = (int*)malloc(3*sizeof(int))) == NULL)
    {
      printf("\nCould not allocate memory for triangle connectivity");
      exit(0);
    }
    /* Read in connectivity */
    /* Indexing should be FORTRAN-like(i.e. numbering starts at 1, instead of 0) */
    fgets(buff,bdim,fp);
    sscanf(buff,"%d %d %d",&tri_conn[t][0],&tri_conn[t][1],&tri_conn[t][2]);
    /* decrement numbers for c indexing */
    --tri_conn[t][0];
    --tri_conn[t][1];
    --tri_conn[t][2];
  }
    
  /* Read number of quadrilateral elements */
  fgets(buff,bdim,fp); // Header text from file
  fgets(buff,bdim,fp); // Line containing number of quad-elements
  sscanf(buff,"%d",&nq);
  printf("\nNumber of quadrilaterals = %d\n",nq);
  if ((quad_conn = (int**)malloc(nq*sizeof(int*))) == NULL)
  {
    printf("\nCould not allocate memory for quadrilateral connectivity");
    exit(0);
  }
  for (q=0; q < nq; q++)
  {
    if ((quad_conn[q] = (int*)malloc(4*sizeof(int))) == NULL)
    {
      printf("\nCould not allocate memory for quadrilateral connectivity");
      exit(0);
    }
    /* Read in connectivity */
    /* Indexing should be FORTRAN-like(i.e. numbering starts at 1, instead of 0) */
    fgets(buff,bdim,fp);
    sscanf(buff,"%d %d %d %d",&quad_conn[q][0],&quad_conn[q][1],
                              &quad_conn[q][2],&quad_conn[q][3]);
    /* decrement numbers for c indexing */
    --quad_conn[q][0];
    --quad_conn[q][1];
    --quad_conn[q][2];
    --quad_conn[q][3];
  }

  fclose(fp);
  
  /* End reading 2D grid file */

  /*
     ELEMENT NODE INDEXING CONVENTIONS
     ---------------------------------
     triangle index notation:
                       2---1        tri_conn[t][0]   t is the number of the triangle
                        \ /         tri_conn[t][1]   t ranges from 0 to nt-1
                         0          tri_conn[t][2]

    Access node numbers as follows:
    n0 = tri_conn[t][0];
    n1 = tri_conn[t][1];
    n2 = tri_conn[t][2];
    Retrieve coordinates as follows:
    x0 = x[n0];
    y0 = y[n0];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];

     quadrilateral notation:
                       3---2        quad_conn[q][0]   q is the number of the quad
                       |   |        quad_conn[q][1]   q ranges from 0 to nq-1
                       0---1        quad_conn[q][2]
                                    quad_conn[q][3]

    Access node numbers as follows:
    n0 = quad_conn[t][0];
    n1 = quad_conn[t][1];
    n2 = quad_conn[t][2];
    n3 = quad_conn[t][3];
    Retrieve coordinates as follows:
    x0 = x[n0];
    y0 = y[n0];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];
    x3 = x[n3];
    y3 = y[n3];

  */

  /* INSERT LOGIC HERE TO PERFORM THE REQUIRED METRIC CALCULATIONS */
  double  x0, x1, x2, x3;
  double  y0, y1, y2, y3;
  double a, amn, avg, amx;
  double ux, uy, vx, vy;
  int n0, n1, n2, n3;

/* === Triangle Area === */

  amn =  1.0e10;
  amx = -1.0e10;
  avg =  0.0;
  if (nt > 0){               // if there are triangles
     for (t=0; t < nt; t++){ // loop over all triangles
	     n0 = tri_conn[t][0];
		  n1 = tri_conn[t][1];
		  n2 = tri_conn[t][2];
		  ux = x[n1]-x[n0];
		  uy = y[n1]-y[n0];
		  vx = x[n2]-x[n0];
		  vy = y[n2]-y[n0];
		  a = 0.5*(ux*vy-uy*vx);
		  amn = MIN(a,amn);
		  amx = MAX(a,amx);
		  avg += a;
	  }
	  printf("\nTriangle minimum area = %17.10e",amn);
	  printf("\nTriangle average area = %17.10e",avg/MAX(1,nt));
	  printf("\nTriangle maximum area = %17.10e\n",amx);
  }

/* === Quad Area === */

  amn =  1.0e10;
  amx = -1.0e10;
  avg =  0.0;
  if (nq > 0){              // if there are quads
	  for (q=0; q < nq; q++) // loop over all quads
	  {
		  n0 = quad_conn[q][0];
		  n1 = quad_conn[q][1];
		  n2 = quad_conn[q][2];
		  n3 = quad_conn[q][3];
		  ux = x[n2]-x[n0];
		  uy = y[n2]-y[n0];
		  vx = x[n3]-x[n1];
		  vy = y[n3]-y[n1];
		  a = 0.5*(ux*vy-uy*vx);
		  amn = MIN(a,amn);
		  amx = MAX(a,amx);
		  avg += a;
	  }
	  printf("\nQuadrilateral minimum area = %17.10e",amn);
	  printf("\nQuadrilateral average area = %17.10e",avg/MAX(1,nq));
	  printf("\nQuadrilateral maximum area = %17.10e\n",amx);
  }

/* === Triangle Angle === */

  double theta, thetamin, thetaavg, thetamax;
  double magu, magv, dotprod;

  thetamin = 1.0e10;
  thetamax = -1.0e10;
  thetaavg = 0.0;
  if (nt > 0){ // if there are tri
	  for (t=0; t < nt; t++){ // loop over all tri
		  for (i = 0; i < 3; i++){ // loop over nodes in tri
          n0 = tri_conn[t][i];
          n1 = tri_conn[t][i+1];
          n2 = tri_conn[t][2*(1-i)];
        if ( i == 2 ){
			 n1 = tri_conn[t][0];
			 n2 = tri_conn[t][1];
		 }
		 ux = x[n1]-x[n0];
		 uy = y[n1]-y[n0];
		 vx = x[n2]-x[n0];
		 vy = y[n2]-y[n0];
	  	 dotprod = ux*vx + uy*vy;
		 magu = sqrt(ux*ux + uy*uy);
		 magv = sqrt(vx*vx + vy*vy);
		 theta = acos(dotprod / (magu*magv));
		 theta = theta * 180.0 / PI;
		 thetamin = MIN(theta,thetamin);
		 thetamax = MAX(theta,thetamax);
		 thetaavg += theta;
		 }
	  }
	  printf("\nTriangle minimum angle = %13.10f degrees",thetamin);
	  printf("\nTriangle average angle = %13.10f degrees",thetaavg/(3*MAX(1,nt)));
	  printf("\nTriangle maximum angle = %13.10f degrees\n",thetamax);
  }

/* === Quad Angle computations === */

  thetamin =  1.0e10;
  thetamax = -1.0e10;
  thetaavg = 0.0;
  if (nq > 0){ // if there are quads
    for (q=0; q < nq; q++){ // loop over all quads
      for (i = 0; i < 4; i++){ // loop over nodes in quad
     	  if (i == 0){
     	  	 n0 = quad_conn[q][0];
          n1 = quad_conn[q][1];
          n2 = quad_conn[q][2];
          n3 = quad_conn[q][3];
     	  }
     	  else if (i == 1){
     	  	 n0 = quad_conn[q][1];
          n1 = quad_conn[q][2];
          n2 = quad_conn[q][3];
          n3 = quad_conn[q][0];
     	  }
     	  else if (i == 2){
     	  	 n0 = quad_conn[q][2];
          n1 = quad_conn[q][3];
          n2 = quad_conn[q][0];
          n3 = quad_conn[q][1];
     	  }
     	  else if (i == 3){ 
     	  	 n0 = quad_conn[q][3];
          n1 = quad_conn[q][0];
          n2 = quad_conn[q][1];
          n3 = quad_conn[q][2];
     	  }
        ux = x[n1]-x[n0];
        uy = y[n1]-y[n0];
        vx = x[n3]-x[n0];
        vy = y[n3]-y[n0];
        dotprod = ux*vx + uy*vy;
        magu = sqrt(ux*ux + uy*uy);
        magv = sqrt(vx*vx + vy*vy);
        theta = acos(dotprod / (magu*magv));
        theta = theta * 180.0 / PI;
        thetamin = MIN(theta,thetamin);
        thetamax = MAX(theta,thetamax);
        thetaavg += theta;
      }
    }
    printf("\nQuadrilateral minimum angle = %14.10f degrees",thetamin);
    printf("\nQuadrilateral average angle = %14.10f degrees",thetaavg/(4*MAX(1,nq)));
    printf("\nQuadrilateral maximum angle = %14.10f degrees\n",thetamax);
  }

  
/* === Triangle AR computation === */

  double ax, ay, la;
  double bx, by, lb;
  double cx, cy, lc;
  double s, AR, ARmin, ARavg, ARmax;

  ARmin = 1.0e10;
  ARmax = -1.0e10;
  ARavg = 0.0;
  if (nt > 0){
    for (t=0; t < nt; t++)
    {
       n0 = tri_conn[t][0];
       n1 = tri_conn[t][1];
       n2 = tri_conn[t][2];
       ax = x[n1]-x[n0];
       ay = y[n1]-y[n0];
       bx = x[n2]-x[n1];
       by = y[n2]-y[n1];
       cx = x[n0]-x[n2];
       cy = y[n0]-y[n2];
       la = sqrt(ax*ax + ay*ay);
       lb = sqrt(bx*bx + by*by);
       lc = sqrt(cx*cx + cy*cy);
       s = 0.5*(la + lb + lc);
       AR = la*lb*lc / (8.0*(s-la)*(s-lb)*(s-lc));
       ARmin = MIN(AR,ARmin);
       ARmax = MAX(AR,ARmax);
       ARavg += AR;
    }
    printf("\nTriangle minimum aspect ratio = %17.10e",ARmin);
    printf("\nTriangle average aspect ratio = %17.10e",ARavg/MAX(1,nt));
    printf("\nTriangle maximum aspect ratio = %17.10e\n",ARmax);
  }

/* === Quad AR computation === */

  double pl1, pl2;
  double mp1x, mp1y;
  double mp2x, mp2y;
  double mp3x, mp3y;
  double mp4x, mp4y;

  ARmin = 1.0e10;
  ARmax = -1.0e10;
  ARavg = 0.0;
  if (nq > 0){ // if there are quads
    for (q=0; q < nq; q++)
    {
       n0 = quad_conn[q][0];
       n1 = quad_conn[q][1];
       n2 = quad_conn[q][2];
       n3 = quad_conn[q][3];
       x0 = x[n0];
       y0 = y[n0];
       x1 = x[n1];
       y1 = y[n1];
       x2 = x[n2];
       y2 = y[n2];
       x3 = x[n3];
       y3 = y[n3];
       
       mp1x = half(x1,x0);
       mp1y = half(y1,y0);

       mp2x = half(x2,x1);
       mp2y = half(y2,y1);

       mp3x = half(x3,x2);
       mp3y = half(y3,y2);

       mp4x = half(x3,x0);
       mp4y = half(y3,y0);

       pl1 = sqrt( (mp3x-mp1x)*(mp3x-mp1x) 
                 + (mp3y-mp1y)*(mp3y-mp1y));

       pl2 = sqrt( (mp2x-mp4x)*(mp2x-mp4x) 
                 + (mp2y-mp4y)*(mp2y-mp4y));

       AR = MAX(pl1,pl2) / MIN(pl1,pl2);
       ARmin = MIN(AR,ARmin);
       ARmax = MAX(AR,ARmax);
       ARavg += AR;
    }
    printf("\nQuadrilateral minimum aspect ratio = %17.10e",ARmin);
    printf("\nQuadrilateral average aspect ratio = %17.10e",ARavg/MAX(1,nq));
    printf("\nQuadrilateral maximum aspect ratio = %17.10e\n",ARmax);
  }

/* === Triangle skewness computation === */

  double skew, skewmin, skewavg, skewmax;
  double lmin, lmax;

  skewmin = 1.0e10;
  skewmax = -1.0e10;
  skewavg = 0.0;
  lmin = 1.0e10;
  lmax = -1.0e10;
  if (nt > 0){
    for (t=0; t < nt; t++)
    {
      n0 = tri_conn[t][0];
      n1 = tri_conn[t][1];
      n2 = tri_conn[t][2];
      ax = x[n1]-x[n0];
      ay = y[n1]-y[n0];
      bx = x[n2]-x[n1];
      by = y[n2]-y[n1];
      cx = x[n0]-x[n2];
      cy = y[n0]-y[n2];
      la = sqrt(ax*ax + ay*ay);
      lb = sqrt(bx*bx + by*by);
      lc = sqrt(cx*cx + cy*cy);

      lmin = MIN(MIN(la,lb),lc);
      lmax = MAX(MAX(la,lb),lc);
      skew = lmax / lmin;
      skewmin = MIN(skew,skewmin);
      skewmax = MAX(skew,skewmax);
      skewavg += skew;

    }
    printf("\nTriangle minimum skewness = %17.10e",skewmin);
    printf("\nTriangle average skewness = %17.10e",skewavg/MAX(1,nt));
    printf("\nTriangle maximum skewness = %17.10e\n",skewmax);
  }

/* ===  Quadrilateral skewness computation === */

  double ldiag1, ldiag2;

  skewmin = 1.0e10;
  skewmax = -1.0e10;
  skewavg = 0.0;
  if (nq > 0){
    for (q=0; q < nq; q++)
    {
      n0 = quad_conn[q][0];
      n1 = quad_conn[q][1];
      n2 = quad_conn[q][2];
      n3 = quad_conn[q][3];
      x0 = x[n0];
      y0 = y[n0];
      x1 = x[n1];
      y1 = y[n1];
      x2 = x[n2];
      y2 = y[n2];
      x3 = x[n3];
      y3 = y[n3];

      ldiag1 = sqrt( (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) );
      ldiag2 = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) );
      skew = MAX(ldiag1,ldiag2) / MIN(ldiag1,ldiag2);
      skewmin = MIN(skew,skewmin);
      skewmax = MAX(skew,skewmax);
      skewavg += skew;
    }
    printf("\nQuadrilateral minimum skewness = %17.10e",skewmin);
    printf("\nQuadrilateral average skewness = %17.10e",skewavg/MAX(1,nq));
    printf("\nQuadrilateral maximum skewness = %17.10e\n",skewmax);
  }

/* === Triangle weighted condition number === */

  double AWI11, AWI12, AWI21, AWI22;
  double WAI11, WAI12, WAI21, WAI22;
  double normAWI, normWAI;
  double r3 = sqrt(3.0);
  double sumsq;
  double WCN, WCNmin, WCNmax, WCNavg;
  double detA;
  
  WCNmin = 1.0e10;
  WCNmax = -1.0e10;
  WCNavg = 0.0;
  if (nt > 0){
    for (t=0; t < nt; t++)
    {
      n0 = tri_conn[t][0];
      n1 = tri_conn[t][1];
      n2 = tri_conn[t][2];
      ux = x[n1]-x[n0];
      uy = y[n1]-y[n0];
      vx = x[n2]-x[n0];
      vy = y[n2]-y[n0];

      AWI11 = vx;
      AWI12 = (2*ux-vx)/r3;
      AWI21 = vy;
      AWI22 = (2*uy-vy)/r3;
      detA  = vx*uy - vy*ux;

      WAI11 = ( uy - 0.5*vy)/detA;
      WAI12 = (-ux + 0.5*vx)/detA;
      WAI21 = -r3*0.5*vy/detA;
      WAI22 =  r3*0.5*vx/detA;

      sumsq = sq(AWI11) + sq(AWI12) + sq(AWI21) + sq(AWI22);
      normAWI = sqrt(sumsq);

      sumsq = sq(WAI11) + sq(WAI12) + sq(WAI21) + sq(WAI22);
      normWAI = sqrt(sumsq);

      WCN = normAWI * normWAI / 2.0;
      WCNmin = MIN(WCN,WCNmin);
      WCNmax = MAX(WCN,WCNmax);
      WCNavg += WCN;
    }
    printf("\nTriangle minimum weighted condition number = %17.10e",WCNmin);
    printf("\nTriangle average weighted condition number = %17.10e",WCNavg/MAX(1,nt));
    printf("\nTriangle maximum weighted condition number = %17.10e\n",WCNmax);
  }

/* === Quad condition number === */

  double CN, CNmin, CNavg, CNmax;
  double A11, A12, A21, A22;
  double AI11, AI12, AI21, AI22;
  double trace;
  double normA, normAI;

  CNmin = 1.0e10;
  CNmax = -1.0e10;
  CNavg = 0.0;
  if (nq > 0){
    for (q=0; q < nq; q++){
      for (i = 0; i < 4; i++){ // loop over nodes in quad
     	  if (i == 0){
     	  	 n0 = quad_conn[q][0];
          n1 = quad_conn[q][1];
          n2 = quad_conn[q][2];
          n3 = quad_conn[q][3];
     	  }
     	  else if (i == 1){
     	  	 n0 = quad_conn[q][1];
          n1 = quad_conn[q][2];
          n2 = quad_conn[q][3];
          n3 = quad_conn[q][0];
     	  }
     	  else if (i == 2){
     	  	 n0 = quad_conn[q][2];
          n1 = quad_conn[q][3];
          n2 = quad_conn[q][0];
          n3 = quad_conn[q][1];
     	  }
     	  else if (i == 3){ 
     	  	 n0 = quad_conn[q][3];
          n1 = quad_conn[q][0];
          n2 = quad_conn[q][1];
          n3 = quad_conn[q][2];
     	  }
        ux = x[n1]-x[n0];
        uy = y[n1]-y[n0];
        vx = x[n3]-x[n0];
        vy = y[n3]-y[n0];

        A11 = vx;	
     	  A12 = ux;
        A21 = vy;
        A22 = uy;
        detA  = vx*uy - vy*ux;

     	  AI11 =  uy/detA;
     	  AI12 = -ux/detA;
     	  AI21 = -vy/detA;
     	  AI22 =  vx/detA;

     	  sumsq = sq(A11) + sq(A12) + sq(A21) + sq(A22);
     	  normA = sqrt(sumsq);
     
     	  sumsq = sq(AI11) + sq(AI12) + sq(AI21) + sq(AI22);
     	  normAI = sqrt(sumsq);
         
     	  CN = normA*normAI/ 2.0;
        CNmin = MIN(CN,CNmin);
        CNmax = MAX(CN,CNmax);
        CNavg += CN;
      }
    }
    printf("\nQuadrilateral minimum condition number = %17.10e",CNmin);
    printf("\nQuadrilateral average condition number = %17.10e",CNavg/(4*MAX(1,nq)));
    printf("\nQuadrilateral maximum condition number = %17.10e\n",CNmax);
  }

/* === Triangle Corner Jacobian === */

  double v1x,v1y,v2x,v2y,vmag;
  int numneg = 0;
  int  ncell = 0;
  int nscell = 0;
  int pscell = 0;
  int  pcell = 0;
  double J,Javg,Jmin, Jmax;

  Jmin = 1.0e10;
  Jmax = -1.0e10;
  if (nt > 0){
    for (t=0; t < nt; t++){
      numneg = 0;
      Javg = 0.0; 
    	for (i = 0; i < 3; i++){ // loop over nodes in tri
        n0 = tri_conn[t][i];
        n1 = tri_conn[t][i+1];
        n2 = tri_conn[t][2*(1-i)];
        if ( i == 2 ){
          n1 = tri_conn[t][0];
          n2 = tri_conn[t][1];
        }
        v1x = x[n1]-x[n0];
        v1y = y[n1]-y[n0];
        v2x = x[n2]-x[n0];
        v2y = y[n2]-y[n0];

        vmag = sqrt(v1x*v1x + v1y*v1y);
        v1x /= vmag;
        v1y /= vmag;
        vmag = sqrt(v2x*v2x + v2y*v2y);
        v2x /= vmag;
        v2y /= vmag;
        J = v1x*v2y - v2x*v1y;
        if(J < 0.0)
           numneg++;
     	  Javg += J;
        Jmin = MIN(Jmin,J);
        Jmax = MAX(Jmax,J);
      }
      Javg /= 3;
      //Determine Jacobian Type
      if (numneg == 3) // negative cell
      	ncell++;
      else if ((numneg < 3) && (Javg < 0.0)) // neg skew
      	nscell++;
      else if ((numneg > 0) && (Javg > 0.0)) // pos skew
      	pscell++;
      else if (numneg == 0) // positive cell
      	pcell++;
      else{
      	printf("Error, Jacobian failure\n");
     	   exit(0);
      }
    }
    printf("\nTriangle stats: number of negative        = %d", ncell);
    printf("\nTriangle stats: number of negative skewed = %d",nscell);
    printf("\nTriangle stats: number of positive skewed = %d",pscell);
    printf("\nTriangle stats: number of positive        = %d", pcell);
    printf("\nMinimum triangle Jacobian = %f", Jmin);
    printf("\nMaximum triangle Jacobian = %f\n", Jmax);
  }


/* === Quad Corner Jacobian === */

  numneg = 0;
  ncell = 0;
  nscell = 0;
  pscell = 0;
  pcell = 0;
  Jmin = 1.0e10;
  Jmax = -1.0e10;
  if (nq > 0){
    for (q=0; q < nq; q++){
      numneg = 0; // count negative Jacobians
      Javg = 0.0; // reset Javg for each tri
      for (i = 0; i < 4; i++){ // loop over nodes in quad
     	  if (i == 0){
     	  	 n0 = quad_conn[q][0];
          n1 = quad_conn[q][1];
          n2 = quad_conn[q][2];
          n3 = quad_conn[q][3];
     	  }
     	  else if (i == 1){
     	  	 n0 = quad_conn[q][1];
          n1 = quad_conn[q][2];
          n2 = quad_conn[q][3];
          n3 = quad_conn[q][0];
     	  }
     	  else if (i == 2){
     	  	 n0 = quad_conn[q][2];
          n1 = quad_conn[q][3];
          n2 = quad_conn[q][0];
          n3 = quad_conn[q][1];
     	  }
     	  else if (i == 3){ 
     	  	 n0 = quad_conn[q][3];
          n1 = quad_conn[q][0];
          n2 = quad_conn[q][1];
          n3 = quad_conn[q][2];
     	  }
        v1x = x[n1]-x[n0];
        v1y = y[n1]-y[n0];
        v2x = x[n3]-x[n0];
        v2y = y[n3]-y[n0];

     	  vmag = sqrt(v1x*v1x + v1y*v1y);
     	  v1x /= vmag;
     	  v1y /= vmag;
     	  vmag = sqrt(v2x*v2x + v2y*v2y);
     	  v2x /= vmag;
     	  v2y /= vmag;
     	  J = v1x*v2y - v2x*v1y;
     	  if(J < 0.0)
          numneg++;
     	  Javg += J;
     	  Jmin = MIN(Jmin,J);
     	  Jmax = MAX(Jmax,J);
      }
      Javg /= 4;
      //Determine Jacobian Type
      if (numneg == 4) // negative cell
      	ncell++;
      else if ((numneg < 4) && (Javg < 0.0)) // neg skew
      	nscell++;
      else if ((numneg > 0) && (Javg > 0.0)) // pos skew
      	pscell++;
      else if (numneg == 0) // positive cell
      	pcell++;
      else{
      	printf("Error, Jacobian failure\n");
     	   exit(0);
      }
    }
    printf("\nQuadrilateral stats: number of negative        = %d", ncell);
    printf("\nQuadrilateral stats: number of negative skewed = %d",nscell);
    printf("\nQuadrilateral stats: number of positive skewed = %d",pscell);
    printf("\nQuadrilateral stats: number of positive        = %d", pcell);
    printf("\nMinimum quadrilateral Jacobian = %f", Jmin);
    printf("\nMaximum quadrilateral Jacobian = %f\n", Jmax);
  }

/* ================= End Metric Calculations Zone =================== */
  /* END OF LOGIC FOR PERFORMING METRIC CALCULATIONS */

  // output Gnuplot file
  char filename[32];
  filename[0]='\0';
  strcat(filename,pArgs[1]);
  strcat(filename,".dat");
  printf("\nFilename = <%s>\n",filename);
  // Open file for write
  if ((fp = fopen(filename,"w")) == NULL)
  {
    printf("\nError opening file <%s>.",filename);
    exit(0);
  }
  for (t=0; t < nt; t++)
  {
    n0 = tri_conn[t][0];
    n1 = tri_conn[t][1];
    n2 = tri_conn[t][2];
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n0],y[n0]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n1],y[n1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n2],y[n2]);
    fprintf(fp,"%19.10e %19.10e 0.0\n\n",x[n0],y[n0]);
  }
  for (q=0; q < nq; q++)
  {
    n0 = quad_conn[q][0];
    n1 = quad_conn[q][1];
    n2 = quad_conn[q][2];
    n3 = quad_conn[q][3];
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n0],y[n0]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n1],y[n1]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n2],y[n2]);
    fprintf(fp,"%19.10e %19.10e 0.0\n",  x[n3],y[n3]);
    fprintf(fp,"%19.10e %19.10e 0.0\n\n",x[n0],y[n0]);
  }
  fclose(fp);

  /* END OF LOGIC FOR PERFORMING METRIC CALCULATIONS */

  /* free up all allocated memory */
  free(x);
  free(y);
  for (t=0; t < nt; t++)
    free(tri_conn[t]);
  free(tri_conn);
  for (q=0; q < nq; q++)
    free(quad_conn[q]);
  free(quad_conn);
}


double half(double p1, double p2)
{
	return (double) ((p1) + (p2)) / 2.0;
}

double sq(double val)
{
   return (double) (val)*(val);
}
