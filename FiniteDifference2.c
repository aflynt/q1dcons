/*
 * File : MacCormack.c 11/29/2012
 * Abstract:
 * Solving the Transient 1-D flow equations using Time Marching and MacCormaks
 * 2nd order scheme.
 */

#define S_FUNCTION_NAME  FiniteDifference2
#define S_FUNCTION_LEVEL 2

/* These define the block's input and output width number of parameters
 * that the block uses and the sample time for Simulink*/
#define OUTPUT_WIDTH 13
#define NPARAM 15
#define WORK 13

/*included files and variables */
#include "simstruc.h"
#include "matrix.h"
#include "mex.h"
#include <stdlib.h> /* This is required to use malloc and free */
#include <string.h>
#include <math.h>

/* Simulink Sfunction Variable List */
#define PARAM_DEF0(S) ssGetSFcnParam(S, 0) /* Tunnel Geom x in m*/
#define PARAM_DEF1(S) ssGetSFcnParam(S, 1) /* Tunnel Geom a in m2*/
#define PARAM_DEF2(S) ssGetSFcnParam(S, 2) /* Source_c terms flow coefficient m2 */
#define PARAM_DEF3(S) ssGetSFcnParam(S, 3) /* Source_m terms loss coefficient */
#define PARAM_DEF4(S) ssGetSFcnParam(S, 4) /* Source_e terms heat transfer area m2*/
#define PARAM_DEF5(S) ssGetSFcnParam(S, 5) /* Section Properties (Length, Alhbc, Arhbc) */
#define PARAM_DEF6(S) ssGetSFcnParam(S, 6) /* Number of Nodes */
#define PARAM_DEF7(S) ssGetSFcnParam(S, 7) /* Node Number of TestSection */
#define PARAM_DEF8(S) ssGetSFcnParam(S, 8) /* ICs (rho,e) */
#define PARAM_DEF9(S) ssGetSFcnParam(S, 9) /* GAS type 0=airl 1=nitrogen */
#define PARAM_DEF10(S) ssGetSFcnParam(S, 10) /* Test Section Grid Range [Ts1 Ts2]*/
#define PARAM_DEF11(S) ssGetSFcnParam(S, 11) /* CDa; Model Drag Coeff (2nd order Poly) */
#define PARAM_DEF12(S) ssGetSFcnParam(S, 12) /* CLa; Model Lift Coeff (2nd order Poly) */
#define PARAM_DEF13(S) ssGetSFcnParam(S, 13) /* Model Frontal Area*/
#define PARAM_DEF14(S) ssGetSFcnParam(S, 14) /* Test Section Loss Adjustment Parameter*/

#define GRID mxGetScalar(PARAM_DEF6(S))

/* Define Function ProtoTypes ********************************************/
void MC_Method(SimStruct *S);
void MC_BC(SimStruct *S);
void MC_SrcTerms(SimStruct *S, real_T q[3], real_T src[3], real_T PLN[2], const int i);
real_T h_convc(const real_T *q, const real_T a, const real_T *gasProp, const real_T Tw, const int EOS);
real_T vis(const real_T T, const int EOS);
/************************************************************************/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{

    ssSetNumSFcnParams(S, NPARAM);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    /* Since this Sfunction does not have any time dependent states these
     * the next two functions define the number of states as zero*/
    ssSetNumContStates(S, 3*GRID);
    ssSetNumDiscStates(S, 0);

    /* Define the number of input ports and their width */
    if (!ssSetNumInputPorts(S, 5)) return;
    ssSetInputPortWidth(S, 0, 3); /* lhbc (rho, mdot, T)*/
    ssSetInputPortRequiredContiguous(S, 0, true);
    ssSetInputPortWidth(S, 1, 3); /* rhbc (rho, P, T)*/
    ssSetInputPortRequiredContiguous(S, 1, true);
    ssSetInputPortWidth(S, 2, 3); /* gasProp (Cv, Cp, R)*/
    ssSetInputPortRequiredContiguous(S, 2, true);
    ssSetInputPortWidth(S, 3, 3); /* Plenum Conditions (P, h, Tw) */
    ssSetInputPortRequiredContiguous(S, 3, true); 
    ssSetInputPortWidth(S, 4, 1); /* Model Angle of Attack AOA (deg) */
    ssSetInputPortRequiredContiguous(S, 4, true); 

    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 0);
    ssSetInputPortDirectFeedThrough(S, 1, 0);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 0);
    ssSetInputPortDirectFeedThrough(S, 4, 0);
    
     /* Define the number of output ports and their width */
    if (!ssSetNumOutputPorts(S, 6)) return;
    ssSetOutputPortWidth(S, 0, GRID); /* rho */
    ssSetOutputPortWidth(S, 1, GRID); /* mdot */
    ssSetOutputPortWidth(S, 2, GRID); /* T static */
    ssSetOutputPortWidth(S, 3, GRID); /* Mach */
    ssSetOutputPortWidth(S, 4, GRID); /* Ps */
    ssSetOutputPortWidth(S, 5, OUTPUT_WIDTH); /* General Outputs */
    
    ssSetNumDWork(S,WORK);
    ssSetDWorkWidth(S, 0, GRID);  /* x   */
    ssSetDWorkWidth(S, 1, GRID);  /* dx+ */
    ssSetDWorkWidth(S, 2, GRID);  /* dx- */
    ssSetDWorkWidth(S, 3, GRID);  /* a   */
    ssSetDWorkWidth(S, 4, GRID);  /* dx  */
    ssSetDWorkWidth(S, 5, GRID);  /* src_c */
    ssSetDWorkWidth(S, 6, GRID);  /* src_m  */
    ssSetDWorkWidth(S, 7, GRID);  /* src_e */
    ssSetDWorkWidth(S, 8, 3);     /* lhnode (rhoA, rU, rhoEA)*/
    ssSetDWorkWidth(S, 9, 3);     /* rhnode (rhoA, rU, rhoEA)*/
    ssSetDWorkWidth(S, 10, 1);    /* sumMdotPL */
    ssSetDWorkWidth(S, 11, 1);    /* sumMdot_hPL  */
    ssSetDWorkWidth(S, 12, GRID); /* dAdx  */
    
    ssSetDWorkDataType(S, 0, SS_DOUBLE);  /* x*/
    ssSetDWorkDataType(S, 1, SS_DOUBLE);  /* dx+*/
    ssSetDWorkDataType(S, 2, SS_DOUBLE);  /* dx-*/
    ssSetDWorkDataType(S, 3, SS_DOUBLE);  /* a*/
    ssSetDWorkDataType(S, 4, SS_DOUBLE);  /* dx  */
    ssSetDWorkDataType(S, 5, SS_DOUBLE);  /* src_c*/
    ssSetDWorkDataType(S, 6, SS_DOUBLE);  /* src_m */
    ssSetDWorkDataType(S, 7, SS_DOUBLE);  /* src_e */
    ssSetDWorkDataType(S, 8, SS_DOUBLE);  /* (rhoA, rU, rhoEA)*/
    ssSetDWorkDataType(S, 9, SS_DOUBLE);  /* (rhoA, rU, rhoEA)*/
    ssSetDWorkDataType(S, 10, SS_DOUBLE);  /* sumMdotPL */
    ssSetDWorkDataType(S, 11, SS_DOUBLE);  /* sumMdot_hPL  */
    ssSetDWorkDataType(S, 12, SS_DOUBLE);  /* dAdx  */

    /* specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S,0);
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    //ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

#define MDL_INITIALIZE_CONDITIONS
void mdlInitializeConditions(SimStruct *S)
{

    /* Continuous States */
    real_T *xC     = ssGetContStates(S);

    /* Parameters */
    const real_T  *secProp= mxGetData(PARAM_DEF5(S));
    const real_T  *N      = mxGetData(PARAM_DEF6(S));
    const real_T  *IC     = mxGetData(PARAM_DEF8(S)); /* IC= [rho,e] */

   /* Work Variables */
   real_T *a           = (real_T*) ssGetDWork(S,3);
   real_T *lhnode      = (real_T*) ssGetDWork(S,8); /* (rhoA, rU, rhoEA)*/
   real_T *rhnode      = (real_T*) ssGetDWork(S,9); /* (rhoA, rU, rhoEA)*/
   real_T *sumMdotPL   = (real_T*) ssGetDWork(S,10);
   real_T *sumMdot_hPL = (real_T*) ssGetDWork(S,11);

   real_T Alhbc = secProp[1];
   real_T Arhbc = secProp[2];
   int_T Nmax = floor(N[0]);
   int_T i;

   for (i=0;i<Nmax;i++){ /*initialize the continuous variables*/
       /*xC[0*Nmax+i] == rA; xC[1*Nmax+i]==rU; xC[2*Nmax+i] == rEA*/
       xC[0*Nmax+i] = IC[0]*a[i];       /* rho A */
       xC[1*Nmax+i] = 0;                /* rho U -quiescent*/
       xC[2*Nmax+i] = IC[0]*a[i]*IC[1]; /* rho A e */

    /* mexPrintf("Made it here Initialize = %f, %i \n",xC[4*Nmax+i],i); //*/

   }
   /* mexPrintf("Made it here Initialize = %f, %i \n",xC[4*Nmax+Nmax-1],Nmax-1); //*/

   sumMdotPL[0] = 0.0;    /*initialize the sum of the mass flow exchange to
                         *plenum*/
   sumMdot_hPL[0] = 0.0;  /*initialize the sum of mdot x enthalpy going into the
                         *plenum*/
   lhnode[0] = IC[0]*Alhbc;
   lhnode[1] = 0;
   lhnode[2] = IC[0]*Alhbc*IC[1];

   rhnode[0] = IC[0]*Arhbc;
   rhnode[1] = 0;
   rhnode[2] = IC[0]*Arhbc*IC[1];

}

#define MDL_START
void mdlStart(SimStruct *S) {
    
   /* Parameters */
   const real_T  *geom_x   = mxGetData(PARAM_DEF0(S));
   const real_T  *geom_a   = mxGetData(PARAM_DEF1(S));
   const real_T  *Source_c = mxGetData(PARAM_DEF2(S));
   const real_T  *Source_m = mxGetData(PARAM_DEF3(S));
   const real_T  *Source_e = mxGetData(PARAM_DEF4(S));
   const real_T  *secProp  = mxGetData(PARAM_DEF5(S));
   const real_T  *N        = mxGetData(PARAM_DEF6(S));
 
  /* Work Variables */
   real_T *x      = (real_T*) ssGetDWork(S,0);
   real_T *dxf     = (real_T*) ssGetDWork(S,1);
   real_T *dxr    = (real_T*) ssGetDWork(S,2);
   real_T *a      = (real_T*) ssGetDWork(S,3);
   real_T *dx     = (real_T*) ssGetDWork(S,4);
   real_T *src_c  = (real_T*) ssGetDWork(S,5);
   real_T *src_m  = (real_T*) ssGetDWork(S,6);
   real_T *src_e  = (real_T*) ssGetDWork(S,7);
   real_T *dAdx   = (real_T*) ssGetDWork(S,12);
      
   real_T xf = secProp[0]; /* Section Length */
   real_T Alhbc = secProp[1]; /* Area at the lhbc */
   real_T Arhbc = secProp[2]; /* Area at the lhbc */
   int_T Nmax = floor(N[0]);
   int i = 0;
   
   for (i=0; i<Nmax; ++i) { /*transfer the input file information to work
                             *variables */
      x[i] = geom_x[i];
      a[i] = geom_a[i];
      src_c[i] = Source_c[i];
      src_m[i] = Source_m[i];
      src_e[i] = Source_e[i];
      

     if (i==Nmax-1){/* Central Difference*/
      dxf[i] = xf - geom_x[i-1];}
     else{
     if (i==Nmax-2){/* Forward Difference*/
      dxf[i] = xf - geom_x[i];}
     else{
      dxf[i] = geom_x[i+2] - geom_x[i];}}
      
     if (i==0){/* Central Difference*/
       dxr[i] = geom_x[i+1];}
     else{
     if (i==1){/* Backward Difference*/
       dxr[i] = geom_x[i];}
     else{
       dxr[i] = geom_x[i] - geom_x[i-2];}}
   
      if (i==0){/* Central Difference*/
       dx[i] = geom_x[i+1];
       dAdx[i] = (geom_a[i+1]-Alhbc)/dx[i];}
       else {
           if (i==Nmax-1){
           dx[i] = xf - geom_x[i-1]; 
           dAdx[i] = (Arhbc-geom_a[i-1])/dx[i];}
           else {
           dx[i] = geom_x[i+1]-geom_x[i-1];
           dAdx[i] = (geom_a[i+1]-geom_a[i-1])/dx[i];}}
      
   /* mexPrintf("Made it here Startup %.6f, %.6f, %i \n",dAdx[i],dx[i],i); //*/
   }

   /* mexPrintf("Made it here Startup %.6f \n",xf - geom_x[Nmax-1]); //*/
 }

#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  /* Additional Comments ======================================================
   */
  static void mdlDerivatives(SimStruct *S)
  {
   /* Insert MacCormak Method Here ***************************************/
   MC_Method(S);
  }
#endif /* MDL_DERIVATIVES

/* Function: mdlOutputs ===================================================
 * Abstract:
 *  
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    int i = 0;
    real_T  rho0,T0,h0,P0,rhoN,PN,TN,hN,mdotN,rhoTS,TTS,MachTS,ToTS,PoTS,PsTS;
    real_T  c,T0o,P0o,TNo,PNo;
    
    /* Continuous States */
    real_T  *xC = ssGetContStates(S); 
    
    /* Output Ports */
    real_T  *y0 = ssGetOutputPortRealSignal(S,0); /*size GRID*/
    real_T  *y1 = ssGetOutputPortRealSignal(S,1); /*size GRID*/
    real_T  *y2 = ssGetOutputPortRealSignal(S,2); /*size GRID*/
    real_T  *y3 = ssGetOutputPortRealSignal(S,3); /*size GRID*/
    real_T  *y4 = ssGetOutputPortRealSignal(S,4); /*size GRID*/
    real_T  *y5 = ssGetOutputPortRealSignal(S,5); /*size OUTPT_WIDTH*/
        
     /* Input Ports */
    const real_T *gasProp = ssGetInputPortRealSignal(S,2); /* Cv, Cp, R */

    /* Work Variables */
    real_T  *x      = (real_T*) ssGetDWork(S,0); /* X-Pos */
    real_T  *a      = (real_T*) ssGetDWork(S,3); /* Area */
    real_T  *mdotPL = (real_T*) ssGetDWork(S,10); /* sumMdotPL */
    real_T  *hPL    = (real_T*) ssGetDWork(S,11);/* sumMdot_hPL  */
    
    /* Parameter Variables */
    const real_T  *secProp = mxGetData(PARAM_DEF5(S)); /* Section Properties (Length, Alhbc, Arhbc) */
    const real_T  *N       = mxGetData(PARAM_DEF6(S)); /* Number of Nodes */
    const real_T  *NTS     = mxGetData(PARAM_DEF7(S)); /* Node Number of TestSection*/
    
    const real_T gam = gasProp[1]/gasProp[0];
    const real_T Rgas = gasProp[2];
    const real_T Cpgas = gasProp[1];
    const real_T Cvgas = gasProp[0];
    const real_T L = secProp[0];
    int_T Nmax = floor(N[0]);
    int_T TSnode = floor(NTS[0]);

/* Begin of Outputs Section *******************************************/
    i = 0;
    for (i=0;i<=Nmax-1;i++){
    y0[i] = xC[0*Nmax+i]/a[i];                                  /* rho */
    y1[i] = xC[1*Nmax+i];                                       /* mdot */
    /* Ts = ((rhoEA-0.5*rhoUA*rhoUA/rhoA)/rhoA)/Cvgas */
    y2[i] = ((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas; /* T Static */
    c = sqrt(gam*Rgas*((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas); 
    y3[i] = xC[1*Nmax+i]/(xC[0*Nmax+i]*c);                                    /* Mach No. */
    y4[i] = pow(1.0+0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/(Cpgas*xC[0*Nmax+i]*xC[0*Nmax+i]*y2[i]),gam/(gam-1.0));  /* Po/P */
    }
                        
    /* Upstream Density */
    i = 0;
    rho0 = xC[0*Nmax+i]/a[i];
    T0 =  ((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas;  
    T0o = T0 + 0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/(Cpgas*xC[0*Nmax+i]*xC[0*Nmax+i]);
    P0 =  rho0*T0*Rgas;
    P0o = P0*pow(T0o/T0,gam/(gam-1.0));
    h0 =  Cpgas*T0o; /*Total Enthalpy*/

    /* Down Stream Density & Temperature */
    i = Nmax-1;
    rhoN = xC[0*Nmax+i]/a[i];
    mdotN = xC[1*Nmax+i];
    TN =  ((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas;  
    TNo = TN + 0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/(Cpgas*xC[0*Nmax+i]*xC[0*Nmax+i]);
    PN =  rhoN*TN*Rgas;
    PNo = PN*pow(TNo/TN,gam/(gam-1.0));
    hN =  Cpgas*TNo; /*Total Enthalpy*/


    /* Test Section Density & Temperature */
    if (TSnode == 0){
        i = 0;}
    else{
        i = TSnode-1;}
    c = sqrt(gam*Rgas*((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas);
    rhoTS = xC[0*Nmax+i]/a[i];
    TTS =  ((xC[2*Nmax+i]
            -0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/xC[0*Nmax+i])/Cvgas;  
    ToTS = TTS + 0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/(Cpgas*xC[0*Nmax+i]*xC[0*Nmax+i]);
    PsTS =  rhoTS*TTS*Rgas;
    PoTS = PsTS*pow(ToTS/TTS,gam/(gam-1.0));
    MachTS = xC[1*Nmax+i]/(xC[0*Nmax+i]*c);
    
    /* Other Desired Outputs Go Here */
    y5[0] = rho0; 
    y5[1] = P0;
    y5[2] = h0;
    y5[3] = rhoN;
    y5[4] = mdotN;
    y5[5] = TNo;
    y5[6] = hN;
    y5[7] = ToTS;
    y5[8] = MachTS;
    y5[9] = mdotPL[0];
    y5[10] = hPL[0];
    y5[11] = P0o-PNo;
    y5[12] = TTS;
  }

/* Function: mdlTerminate =================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
}

/*  Other Subfunctions ***************************************************/

void MC_Method(SimStruct *S) 
{
    int i;    
    
    /* Continuous States and Derivative*/
    real_T *dxC = ssGetdX(S);
    real_T *xC = ssGetContStates(S);

    /* Input Variables */
    const real_T *lhbc =  ssGetInputPortRealSignal(S,0); /*(rho, mdot, T)*/
    const real_T *rhbc = ssGetInputPortRealSignal(S,1);  /*(rho, P, T)*/
    const real_T *gasProp = ssGetInputPortRealSignal(S,2); /* gasProp (Cv, Cp, R)*/
    const real_T gam = gasProp[1]/gasProp[0];
    
    /* Work Variables */
    real_T  *x       = (real_T*) ssGetDWork(S,0); /* X-Pos */
    real_T  *dxf     = (real_T*) ssGetDWork(S,1); /* dx */
    real_T  *dxr     = (real_T*) ssGetDWork(S,2); /* dx */
    real_T  *a       = (real_T*) ssGetDWork(S,3); /* a */
    real_T *lhnode   = (real_T*) ssGetDWork(S,8); /* (rhoA, rU, rhoEA)*/
    real_T *rhnode   = (real_T*) ssGetDWork(S,9); /* (rhoA, rU, rhoEA)*/
    real_T *sumMdotPL   = (real_T*) ssGetDWork(S,10); /* sumMdotPL */
    real_T *sumMdot_hPL = (real_T*) ssGetDWork(S,11); /* sumMdot_hPL */
    real_T *dAdx     = (real_T*) ssGetDWork(S,12); /* dAdx */
    
    /* Parameters */
    const real_T  *secProp  = mxGetData(PARAM_DEF5(S));
    const real_T  *N        = mxGetData(PARAM_DEF6(S));
    
    real_T Alhbc = secProp[1]; /*Area of lhbc */
    real_T Arhbc = secProp[2]; /*Area of rhbc */
    
    int_T Nmax = floor(N[0]);
    real_T src[3]={0};
    real_T PLN[2] ={0};
    real_T q[3] = {0};
    real_T Value;
        
    MC_BC(S);  /* set the boundary conditions */
    
   /* MacCormak Method */
    for (i=0;i<Nmax;i++){

   /* Mixed Difference *************************************/ 
    q[0] = xC[0*Nmax+i]; 
    q[1] = xC[1*Nmax+i];
    q[2] = xC[2*Nmax+i];
    
     MC_SrcTerms(S,q,src,PLN,i); /* set the source terms */
   /* Continuity Equation************************************************/
    if (i==0){ /*Central Difference */
     dxC[0*Nmax+i] = (src[0]-(xC[1*Nmax+i+1]-lhnode[1]))/dxr[i];
     }
    else{/* Backward Difference */
     if (i==1){/*Second Order Difference */
     dxC[0*Nmax+i] = (src[0]-(3.0*xC[1*Nmax+i]-4.0*xC[1*Nmax+i-1]+lhnode[1]))/dxr[i];}
     else{
     dxC[0*Nmax+i] = (src[0]-(3.0*xC[1*Nmax+i]-4.0*xC[1*Nmax+i-1]+xC[1*Nmax+i-2]))/dxr[i];}}
        
     /*mexPrintf("Made it here MC_Method = %f, %i \n",rhnode[1],i); //*/
     
   /* Momentum Equation**************************************************/
   if (i==Nmax-1){/*Central Difference */
     dxC[1*Nmax+i] = (src[1]
     /*   d(rUA)/dt = -d(rUA*rUA/rA)/dx - Ai d((gam-1)*(reA-0.5*rUA*rUA/rA)/A)/dx */
                      -((rhnode[1]*rhnode[1]/rhnode[0]+(gam-1.0)*a[i]*(rhnode[2]-0.5*rhnode[1]*rhnode[1]/rhnode[0])/Arhbc)
                       -(xC[1*Nmax+i-1]*xC[1*Nmax+i-1]/xC[0*Nmax+i-1]+(gam-1.0)*a[i]*(xC[2*Nmax+i-1]-0.5*xC[1*Nmax+i-1]*xC[1*Nmax+i-1]/xC[0*Nmax+i-1])/a[i-1])
                       )
                       )/dxf[i];}
     
   /* (src[1]-((xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])-(lhnode[1]*lhnode[1]/lhnode[0]))
                            -a[i]*(gam-1.0)*(((xC[2*Nmax+i]-0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/a[i])
                                            -((lhnode[2]-0.5*lhnode[1]*lhnode[1]/lhnode[0])/Alhbc)))/dxr[i]*/
     
    else{/* Forward Difference */
      if (i==Nmax-2){/*Second Order Difference */
      dxC[1*Nmax+i] = (src[1]
     /*   d(rUA)/dt = -d(rUA*rUA/rA)/dx - Ai d((gam-1)*(reA-0.5*rUA*rUA/rA)/A)/dx */
                      -(-1.0*(rhnode[1]*rhnode[1]/rhnode[0]+(gam-1.0)*a[i]*(rhnode[2]-0.5*rhnode[1]*rhnode[1]/rhnode[0])/Arhbc)
                        +4.0*(xC[1*Nmax+i+1]*xC[1*Nmax+i+1]/xC[0*Nmax+i+1]+(gam-1.0)*a[i]*(xC[2*Nmax+i+1]-0.5*xC[1*Nmax+i+1]*xC[1*Nmax+i+1]/xC[0*Nmax+i+1])/a[i+1])
                        -3.0*(xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i]+(gam-1.0)*a[i]*(xC[2*Nmax+i]-0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/a[i])
                       )
                       )/dxf[i];}
         else{
      dxC[1*Nmax+i] = (src[1]
     /*   d(rUA)/dt = -d(rUA*rUA/rA)/dx - Ai d((gam-1)*(reA-0.5*rUA*rUA/rA)/A)/dx */
                      -(-1.0*(xC[1*Nmax+i+2]*xC[1*Nmax+i+2]/xC[0*Nmax+i+2]+(gam-1.0)*a[i]*(xC[2*Nmax+i+2]-0.5*xC[1*Nmax+i+2]*xC[1*Nmax+i+2]/xC[0*Nmax+i+2])/a[i+2])
                        +4.0*(xC[1*Nmax+i+1]*xC[1*Nmax+i+1]/xC[0*Nmax+i+1]+(gam-1.0)*a[i]*(xC[2*Nmax+i+1]-0.5*xC[1*Nmax+i+1]*xC[1*Nmax+i+1]/xC[0*Nmax+i+1])/a[i+1])
                        -3.0*(xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i]+(gam-1.0)*a[i]*(xC[2*Nmax+i]-0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])/a[i])
                       )
                       )/dxf[i];}}
    
     /*mexPrintf("Made it here MC_Method = %f, %f \n",dxC[1*Nmax+0],dxC[1*Nmax+1]); //*/
   
   /* Energy Equation****************************************************/
   if (i==0){ /*Central Difference */
     dxC[2*Nmax+i] = (src[2]-(((xC[1*Nmax+i+1]/xC[0*Nmax+i+1])*(xC[2*Nmax+i+1]+(gam-1.0)*(xC[2*Nmax+i+1]-0.5*xC[1*Nmax+i+1]*xC[1*Nmax+i+1]/xC[0*Nmax+i+1])))
                             -((lhnode[1]/lhnode[0])*(lhnode[2]+(gam-1.0)*(lhnode[2]-0.5*lhnode[1]*lhnode[1]/lhnode[0])))))/dxr[i];}
    else{ /* Backward Difference */
     if (i==1){ /*Second Order Difference */
     dxC[2*Nmax+i] = (src[2]-(3.0*((xC[1*Nmax+i]/xC[0*Nmax+i])*(xC[2*Nmax+i]+(gam-1.0)*(xC[2*Nmax+i]-0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])))
                             -4.0*((xC[1*Nmax+i-1]/xC[0*Nmax+i-1])*(xC[2*Nmax+i-1]+(gam-1.0)*(xC[2*Nmax+i-1]-0.5*xC[1*Nmax+i-1]*xC[1*Nmax+i-1]/xC[0*Nmax+i-1])))
                             +1.0*((lhnode[1]/lhnode[0])*(lhnode[2]+(gam-1.0)*(lhnode[2]-0.5*lhnode[1]*lhnode[1]/lhnode[0])))))/dxr[i];}
     else{
     dxC[2*Nmax+i] = (src[2]-(3.0*((xC[1*Nmax+i]/xC[0*Nmax+i])*(xC[2*Nmax+i]+(gam-1.0)*(xC[2*Nmax+i]-0.5*xC[1*Nmax+i]*xC[1*Nmax+i]/xC[0*Nmax+i])))
                             -4.0*((xC[1*Nmax+i-1]/xC[0*Nmax+i-1])*(xC[2*Nmax+i-1]+(gam-1.0)*(xC[2*Nmax+i-1]-0.5*xC[1*Nmax+i-1]*xC[1*Nmax+i-1]/xC[0*Nmax+i-1])))
                             +1.0*((xC[1*Nmax+i-2]/xC[0*Nmax+i-2])*(xC[2*Nmax+i-2]+(gam-1.0)*(xC[2*Nmax+i-2]-0.5*xC[1*Nmax+i-2]*xC[1*Nmax+i-2]/xC[0*Nmax+i-2])))))/dxr[i];}}
  
   
    //mexPrintf("Made it here MC_Method = %f, \t %f, \t %f \t %i \n",PLN[0],PLN[1],dxC[0*Nmax+i],i);
     
   } /* END OF LOOP******************************************************/

   sumMdotPL[0]   =  PLN[0];
   sumMdot_hPL[0] =  PLN[1];
    
}

/************************************************************************/
void MC_BC(SimStruct *S) /* Set the boundary conditions */
{
   real_T *xC = ssGetContStates(S);
   
   /* Work Variables located in the SimStruct S */
   real_T *x      = (real_T*) ssGetDWork(S,0);
   real_T *a      = (real_T*) ssGetDWork(S,3);
   real_T *lhnode  = (real_T*) ssGetDWork(S,8); /* (rhoA, mdot, rhoEA)*/
   real_T *rhnode  = (real_T*) ssGetDWork(S,9); /* (rhoA, mdot, rhoEA)*/
   
   /* Parameters */
    const real_T  *secProp  = mxGetData(PARAM_DEF5(S));
    const real_T  *N        = mxGetData(PARAM_DEF6(S));
   
   /* Input Variables */
   /* port 0 upstream (lhbc) */
    const real_T *lhbc = ssGetInputPortRealSignal(S,0); /*(rho, mdot, T)*/
   /* port 1 dnstream (rlhbc) */
    const real_T *rhbc = ssGetInputPortRealSignal(S,1);/*(rho, P, T)*/
   /* port 2 gas properties) */
    const real_T *gasProp = ssGetInputPortRealSignal(S,2); /*Cv,Cp,R*/
    const real_T gam = gasProp[1]/gasProp[0];
   
   real_T mach,u,p,pb;
   real_T Alhbc = secProp[1];
   real_T Arhbc = secProp[2];
   int_T Nmax = floor(N[0]);
   int_T status;
   
   /* lhbc **************************************************************/ 
   /*the lhbc is applied on the corrector step use predictor values
    *therefore use xC[0*Nmax+0] = rA; xC[1*Nmax+0] = rUA; xC[2*Nmax+0] = rEA*/   
    
      mach=fabs((lhbc[1]/(lhnode[0]))/(pow(gam*gasProp[2]*lhbc[2],0.5)));
      
      if (mach >= 1.0) {
         lhnode[0] = lhbc[0]*Alhbc;
         lhnode[1] = lhbc[1];
         u = lhbc[1] / lhnode[0];
         lhnode[2] = lhnode[0]*(gasProp[0]*lhbc[2]+0.5*u*u);} 
      else { 
         /*Extrapolate rho and fix (rUA) and e*/
         lhnode[0] = lhbc[0]*Alhbc; //xC[0*Nmax+0]*Alhbc/a[0];
         lhnode[1] = lhbc[1];
         u = lhbc[1] / lhnode[0];
         lhnode[2] = lhnode[0]*(gasProp[0]*lhbc[2]+0.5*u*u);}
      
   /* rhbc **************************************************************/ 
   /*the rhbc is applied on the predictor step use corrector values
    *therefore use xC[3*Nmax+Nmax-1] = rA; xC[4*Nmax+Nmax-1] = rU; xC[5*Nmax+Nmax-1] = rEA */
     if ( xC[1*Nmax+Nmax-1] > 0 ) { /*Outgoing flow*/
         rhnode[0]=xC[0*Nmax+Nmax-1]*Arhbc/a[Nmax-1];
         rhnode[1]=xC[1*Nmax+Nmax-1];
         u = rhnode[1]/rhnode[0];
         p = (gam-1)*(xC[2*Nmax+Nmax-1]-0.5*xC[1*Nmax+Nmax-1]*xC[1*Nmax+Nmax-1]/xC[0*Nmax+Nmax-1])/a[Nmax-1];
         pb = rhbc[1];
         rhnode[2]=Arhbc*(xC[2*Nmax+Nmax-1]/a[Nmax-1]+(pb-p)/(gam-1));}
      else { /*Incoming flow; condition*/
         rhnode[0]=rhbc[0]*Arhbc;
         rhnode[1]=xC[1*Nmax+Nmax-1]; 
         u = rhnode[1]/rhnode[0];
         rhnode[2]=rhbc[0]*gasProp[0]*rhbc[2]*Arhbc + 0.5*rhbc[0]*Arhbc*u*u;}
                  
    /*  Non-Reflective BC */
   /*  rhnode[0]=xC[0*Nmax+Nmax-1]*Arhbc/a[Nmax-1];
      rhnode[1]=xC[1*Nmax+Nmax-1];
      rhnode[2]=xC[2*Nmax+Nmax-1]*Arhbc/a[Nmax-1];/**/
      
     /* mexPrintf("Made it here MC_BC = %f, %f, %i \n",lhbc[1],lhnode[1],0); //*/
     }

/*************************************************************************/
void MC_SrcTerms(SimStruct *S, real_T q[3], real_T src[3], real_T PLN[2], const int i)
{
  
   /* Work Variables*/
   real_T *x    = (real_T*) ssGetDWork(S,0);
   real_T *a    = (real_T*) ssGetDWork(S,3);
   real_T *src_c = (real_T*) ssGetDWork(S,5); /* flow coefficients for mass flow out */
   real_T *src_m = (real_T*) ssGetDWork(S,6); /* loss coefficients for internal momentum loss */
   real_T *src_e = (real_T*) ssGetDWork(S,7); /* conduction heat transfer area out to plenum */
   
   /* Parameters Variables */
   const real_T *N        = mxGetData(PARAM_DEF6(S)); /* Node Number of TestSection */
   const real_T *TS       = mxGetData(PARAM_DEF7(S)); /* Node Number of TestSection */
   const real_T *GAS      = mxGetData(PARAM_DEF9(S)); /* GAS type 0=airl 1=nitrogen */
   const real_T *TSection = mxGetData(PARAM_DEF10(S)); /* Test Section Grid Range */
   const real_T *CDa      = mxGetData(PARAM_DEF11(S)); /* CDa; Model Drag Coeff (2nd order Poly) */
   const real_T *CLa      = mxGetData(PARAM_DEF12(S)); /* CLa; Model Lift Coeff (2nd order Poly) */
   const real_T *ModelArea = mxGetData(PARAM_DEF13(S));/* Model Frontal Area*/
   const real_T *LossPara = mxGetData(PARAM_DEF14(S)); /* Test Section Loss Adjustment Parameter*/
   
   /* Input Variables*/
   const real_T *gasProp = ssGetInputPortRealSignal(S,2);/* (Cv, Cp, R)*/
   const real_T *Plenum = ssGetInputPortRealSignal(S,3); /* (P, h, Tw)*/
   const real_T *AOA = ssGetInputPortRealSignal(S,4); /* Deg */  
      
   const real_T gam = gasProp[1]/gasProp[0];
   const real_T Cvgas = gasProp[0];
   const real_T Cpgas = gasProp[1];
   const real_T Rgas = gasProp[2];
   const int_T TS1 = floor(TSection[0]);
   const int_T TS2 = floor(TSection[1]);
   const int_T TSnode = floor(TS[0]);
   const int_T Nmax = floor(N[0]);
   const int_T EOS = floor(GAS[0]);
     
   const real_T  pi=3.141592653589793238462643383279502884197169399375;
    
   /* Initialize the source terms */
   src[0] = 0; /* mdot ... kg/s */
   src[1] = 0; /* Force/Area ... N/m2 */
   src[2] = 0; /* Power ... J/s */  
   
   { /* for mass flow source ********************************************/
         real_T h_src, sign, C;
         real_T uj = q[1]/q[0]; /* calc velocity */
         real_T tem= (q[2]-0.5*q[1]*q[1]/q[0])/q[0]/Cvgas;/* calc temperature */
         real_T h=gasProp[1]*tem; /* calc enthalpy */
         real_T p = q[0]*Rgas*tem/(a[i]); /* calc pressure */
         
          /* Flow coefficient for mass flow out */
         if (i==0 || i==Nmax-1){
         C=src_c[i];}
         else{
         C = ((x[i]-x[i-1])*((src_c[i-1]+src_c[i]))/(x[i+1]-x[i-1])
                     +(x[i+1]-x[i])*((src_c[i]+ src_c[i+1]))/(x[i+1]-x[i-1]));} /* C => src_c has units of m2 */        

         sign = Plenum[0]-p < 0 ? -1 : 1; /* determine dirction of mass flow */
     
         /* Continuity Portion*/
         src[0] += C*sign*sqrt(fabs(Plenum[0]-p)*q[0]/a[i]); /* kg/s */
         PLN[0] += src[0];
         
         /* Momentum Portion*/
         if (sign < 0) {
           src[1] += (src[0] * uj);} /* N */
         else{
           src[1] += 0;}  
         
         /* Energy Portion */
         h_src = sign < 0 ? h : Plenum[1];
         src[2] += src[0] * (h_src); /* J/s */
         PLN[1] += src[2];
          
      }
   
   /* Loss Coefficents *****************************************/
     { 
       real_T uj = q[1]/q[0]; /* calc velocity */
       real_T tem=(q[2]-0.5*q[1]*q[1]/q[0])/q[0]/Cvgas;/* calc temperature */
       real_T c=sqrt(gam*tem*Rgas);
       real_T Dia = sqrt(a[i]*4.0/pi);
       real_T Tt = tem+0.5*q[1]*q[1]/q[0]/q[0]/Cpgas;
       real_T Kto, CRe, Ktbase,Loss,Red,mu,Pt,p, Mach;
       real_T Ktmodel,CD,CL;
       
       Mach = fabs(q[1]/q[0]/c);
       
       p = q[0]*Rgas*tem/a[i];
       Pt = p*pow(Tt/tem,gam/(gam-1.0));
       mu = vis(tem,EOS);       
       Red = pow((1.0+0.5*(gam-1.0)*Mach*Mach),(-gam/(gam-1.0))+0.5)*Pt*Dia*sqrt(gam/(Rgas*Tt))/mu;
       
       /*Use the value in the log file for other locations than those that follow */
       if (i==0 || i==Nmax-1){
       Loss = src_m[i];}
       else{
        Loss = ((x[i]-x[i-1])*((src_m[i-1]+src_m[i]))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*((src_m[i]+ src_m[i+1]))/(x[i+1]-x[i-1]));}
       
      /* Model Loss **********************************************/
       if (i==TSnode-1){
       CD = (CDa[0] + CDa[1]*Mach + CDa[2]*Mach*Mach);
       CL = (CLa[0] + CLa[1]*Mach + CLa[2]*Mach*Mach);
       Ktmodel = (ModelArea[0]/a[TSnode])*(CD*cos(AOA[0]*pi/180.0)+CL*sin(AOA[0]*pi/180.0));
       Loss = ((x[i]-x[i-1])*((src_m[i-1]+src_m[i]))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*((src_m[i]+ (Ktmodel+src_m[i+1])))/(x[i+1]-x[i-1]));
       }
       
       if (i==TSnode){
       CD = (CDa[0] + CDa[1]*Mach + CDa[2]*Mach*Mach);
       CL = (CLa[0] + CLa[1]*Mach + CLa[2]*Mach*Mach);
       Ktmodel = (ModelArea[0]/a[TSnode])*(CD*cos(AOA[0]*pi/180.0)+CL*sin(AOA[0]*pi/180.0));
       Loss = ((x[i]-x[i-1])*((src_m[i-1]+(Ktmodel + src_m[i])))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*(((Ktmodel+src_m[i])+ src_m[i+1]))/(x[i+1]-x[i-1]));
        }
       
       if (i==TSnode+1){
       CD = (CDa[0] + CDa[1]*Mach + CDa[2]*Mach*Mach);
       CL = (CLa[0] + CLa[1]*Mach + CLa[2]*Mach*Mach);
       Ktmodel = (ModelArea[0]/a[TSnode])*(CD*cos(AOA[0]*pi/180.0)+CL*sin(AOA[0]*pi/180.0));
       Loss = ((x[i]-x[i-1])*(((Ktmodel+src_m[i-1])+src_m[i]))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*((src_m[i]+ src_m[i+1]))/(x[i+1]-x[i-1]));
       }
       
       /*****************************************************************/

       /* Empty Tunnel Loss ********************************************/
       if (i==TS2-3){
         Kto = 0.0911*Mach*Mach*Mach+0.016*Mach*Mach-0.0088*Mach+0.051+LossPara[0];
         //Ktbase (test leg corrected for Red, empty tunnel
         //Reynolds number effect
         if (Mach <= 1){
           CRe = 0.023;}
         else{
            CRe = 0.0517*Mach - 0.0263;}
         
        Ktbase = (Kto - CRe * log(Red / 1.0e7) / log(10.0));
                //Test Leg Loss
        Loss = ((x[i]-x[i-1])*((src_m[i-1]+src_m[i]))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*((src_m[i]+ (Ktbase+src_m[i+1])))/(x[i+1]-x[i-1]));
       }
       
       if (i==TS2-2){
        Kto = 0.0911*Mach*Mach*Mach+0.016*Mach*Mach-0.0088*Mach+0.051+LossPara[0];
         //Ktbase (test leg corrected for Red, empty tunnel
         //Reynolds number effect
         if (Mach <= 1){
           CRe = 0.023;}
         else{
            CRe = 0.0517*Mach - 0.0263;}
         
        Ktbase = (Kto - CRe * log(Red / 1.0e7) / log(10.0));
                //Test Leg Loss
        
        Loss = ((x[i]-x[i-1])*((src_m[i-1]+(Ktbase + src_m[i])))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*(((Ktbase+src_m[i])+ src_m[i+1]))/(x[i+1]-x[i-1]));
       }

       if (i==TS2-1){
        Kto = 0.0911*Mach*Mach*Mach+0.016*Mach*Mach-0.0088*Mach+0.051+LossPara[0];
         //Ktbase (test leg corrected for Red, empty tunnel
         //Reynolds number effect
         if (Mach <= 1){
           CRe = 0.023;}
         else{
            CRe = 0.0517*Mach - 0.0263;}
         
        Ktbase = (Kto - CRe * log(Red / 1.0e7) / log(10.0));
                //Test Leg Loss
        Loss = ((x[i]-x[i-1])*(((Ktbase+src_m[i-1])+src_m[i]))/(x[i+1]-x[i-1])
              +(x[i+1]-x[i])*((src_m[i]+ src_m[i+1]))/(x[i+1]-x[i-1]));
       }

     /*******************************************************************/
       
       /* Kt Loss = - Kt
       // The above loss values are calculated outside of this module 
       // and their values are put into the matrix src_m work vector.
       // It is assumed that frction is constant even though it is a function
       // of Reynolds number.  Plus it is small compared to Kt losses. */
       src[1] += -((Loss)/2.0)*q[1]*fabs(q[1]/q[0]); /* N */
       
       //mexPrintf("Value is %.5f node %i\n",src_m[i],i);
       }
      
   { /*For energy source*/
     
       real_T tem_avg=(q[2]-0.5*q[1]*q[1]/q[0])/q[0]/Cvgas;/* calc temperature */
       real_T tem_inf = Plenum[2]; /* Temperature at Plenum Wall */
       real_T h1, h, h2,ht;
        // The following calculates the individual overall heat transfer coefficient at each node
        // the value src_e is 1 or 0 depending on whether HT occurs
       if(i==0 || i==Nmax-1){
       ht = src_e[i]*h_convc(q,a[i],gasProp,tem_inf,EOS);  
       }
       else {   
        h1 = src_e[i-1]*h_convc(q,a[i-1],gasProp,tem_inf,EOS);
        h  = src_e[i]*h_convc(q,a[i],gasProp,tem_inf,EOS);
        h2 = src_e[i+1]*h_convc(q,a[i+1],gasProp,tem_inf,EOS);
        
        // Calculates the heat transfer coefficient applied at the node
        ht = (x[i]-x[i-1])*(h1+h)/(x[i+1]-x[i-1]) + (x[i+1]-x[i])*(h+h2)/(x[i+1]-x[i-1]);} // J/s/K
       
        /* heat transfer 
        // ht = h*Asurface -> src_e is surface area */
        src[2] += (tem_inf - tem_avg)*ht; /* J/s */
     } 
    
   /* mexPrintf("Value is %.5f \n",src[2]); //*/
   
  } /******************************* END OF FUNCTION ******************/


real_T h_convc(const real_T *q, const real_T a, const real_T *gasProp, const real_T Tw, const int EOS) 
{ /*q[0]=rho*A, q[1]=rho*u*A; q[2]=rho*e*A; a: area; gasProp[0] = Cv; gasProp[1] = Cp; gasProp[2] = R */
   
   const real_T  Tref=273;        /* Kelvin */
   const real_T  pi=3.141592653589793238462643383279502884197169399375;
    
   /*Air*/
   const real_T  vis_oA=1.716e-5;  // kg/(m s)
   const real_T  CA = 111;         // Kelvin
   const real_T  k_oA =0.0241;     // W/(m K)
   const real_T  DA = 194;         // Kelvin
   const real_T  MW_Air = 28.9644;   // Air Molecular Weight
 
   /*Nitrogen*/
   const real_T  vis_oN2 = 1.6630e-5; //(kg/(m s))
   const real_T  CN2 = 107;
   const real_T  k_oN2 =0.0242;     // W/(m K)
   const real_T  DN2 = 150;         // Kelvin
   const real_T  MW_N2 = 28.0134;   // Air Molecular Weight
 
   const real_T gam = gasProp[1]/gasProp[0];
   const real_T Cvgas = gasProp[0];
   const real_T Cpgas = gasProp[1];
   const real_T Rgas = gasProp[2];
   
   const real_T uj = q[1]/q[0]; /* calc velocity */
   const real_T Ta= (q[2]-0.5*q[1]*q[1]/q[0])/q[0]/Cvgas;/* calc temperature */
   
   real_T vis_o, C, k_o, D, vis_gas, vis_w, k_gas, d, Pr, Re;

    if (EOS == 1){
     vis_o = vis_oN2;
     C = CN2;
     k_o = k_oN2;
     D = DN2;}
    else{
     vis_o = vis_oA;
     C = CA;
     k_o = k_oA;
     D = DA;}
   
    /* Viscosity - Sutherland Formula (kg /(m s)) */
    vis_gas = (vis_o)*pow((Ta/Tref),1.5)*((Tref + C)/(Ta + C));
    vis_w = (vis_o)*pow((Tw/Tref),1.5)*((Tref + C)/(Tw + C));
    
    /* Thermal Conductivity - Sutherland Formula (W/(m K))*/
    k_gas =  (k_o)*pow((Ta/Tref),1.5)*((Tref + D)/(Ta + D));
    
   d=sqrt(4.0*a/pi);
   Pr = Cpgas*vis_gas/k_gas;
   Re = fabs(q[1]*d/vis_gas);
   /*Sieder and Tate equation*/
   return 0.027*(k_gas/d)*pow(Re,0.8)*pow(Pr,1.0/3.0)*pow(vis_gas/vis_w,0.14); 
   
}

real_T vis(const real_T T, const int EOS) 
{
  
 /* Define various constants for use in Sutherland Viscosity and Thermal
  * Conductivity Formulas
  */
  const real_T  Tref=273;        // Kelvin

  //Air
  const real_T  vis_oA=1.716e-5;  // kg/(m s)
  const real_T  CA = 111;         // Kelvin
 
 //Nitrogen
  const real_T  vis_oN2 = 1.6630e-5; //(kg/(m s))
  const real_T  CN2 = 107;
  
 real_T vis_o, k_o, C,D;
  
  if (EOS == 1){
     vis_o = vis_oN2;
     C = CN2;
  }
  else{
     vis_o = vis_oA;
     C = CA;
  } 
  
 /* Viscosity - Sutherland Formula (kg /(m s)) */
 return (vis_o)*pow((T/Tref),1.5)*((Tref + C)/(T + C));
}

/***********************************************************************/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
