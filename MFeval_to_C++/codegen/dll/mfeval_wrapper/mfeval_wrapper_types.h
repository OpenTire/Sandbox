/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mfeval_wrapper_types.h
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

#ifndef MFEVAL_WRAPPER_TYPES_H
#define MFEVAL_WRAPPER_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_b_struct_T
#define typedef_b_struct_T

typedef struct {
  boolean_T useLimitsCheck;
  boolean_T useAlphaStar;
  boolean_T useTurnSlip;
  boolean_T isLowSpeed;
  boolean_T isLowSpeedAlpha;
  double userDynamics;
} b_struct_T;

#endif                                 /*typedef_b_struct_T*/

#ifndef typedef_f_struct_T
#define typedef_f_struct_T

typedef struct {
  double Fx;
  double Fy;
  double Fz;
  double Mx;
  double My;
  double Mz;
} f_struct_T;

#endif                                 /*typedef_f_struct_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  char FILE_TYPE[5];
  double FILE_VERSION;
  char FILE_FORMAT[7];
  char LENGTH[7];
  char FORCE[8];
  char ANGLE[9];
  char MASS[4];
  char TIME[8];
  double FITTYP;
  char TYRESIDE[6];
  double LONGVL;
  double VXLOW;
  double ROAD_INCREMENT;
  double ROAD_DIRECTION;
  char PROPERTY_FILE_FORMAT[6];
  double USER_SUB_ID;
  double N_TIRE_STATES;
  double USE_MODE;
  double HMAX_LOCAL;
  double TIME_SWITCH_INTEG;
  double UNLOADED_RADIUS;
  double WIDTH;
  double ASPECT_RATIO;
  double RIM_RADIUS;
  double RIM_WIDTH;
  double INFLPRES;
  double NOMPRES;
  double MASS1;
  double IXX;
  double IYY;
  double BELT_MASS;
  double BELT_IXX;
  double BELT_IYY;
  double GRAVITY;
  double FNOMIN;
  double VERTICAL_STIFFNESS;
  double VERTICAL_DAMPING;
  double MC_CONTOUR_A;
  double MC_CONTOUR_B;
  double BREFF;
  double DREFF;
  double FREFF;
  double Q_RE0;
  double Q_V1;
  double Q_V2;
  double Q_FZ2;
  double Q_FCX;
  double Q_FCY;
  double Q_CAM;
  double PFZ1;
  double BOTTOM_OFFST;
  double BOTTOM_STIFF;
  double LONGITUDINAL_STIFFNESS;
  double LATERAL_STIFFNESS;
  double YAW_STIFFNESS;
  double FREQ_LONG;
  double FREQ_LAT;
  double FREQ_YAW;
  double FREQ_WINDUP;
  double DAMP_LONG;
  double DAMP_LAT;
  double DAMP_YAW;
  double DAMP_WINDUP;
  double DAMP_RESIDUAL;
  double DAMP_VLOW;
  double Q_BVX;
  double Q_BVT;
  double PCFX1;
  double PCFX2;
  double PCFX3;
  double PCFY1;
  double PCFY2;
  double PCFY3;
  double PCMZ1;
  double Q_RA1;
  double Q_RA2;
  double Q_RB1;
  double Q_RB2;
  double ELLIPS_SHIFT;
  double ELLIPS_LENGTH;
  double ELLIPS_HEIGHT;
  double ELLIPS_ORDER;
  double ELLIPS_MAX_STEP;
  double ELLIPS_NWIDTH;
  double ELLIPS_NLENGTH;
  double PRESMIN;
  double PRESMAX;
  double FZMIN;
  double FZMAX;
  double KPUMIN;
  double KPUMAX;
  double ALPMIN;
  double ALPMAX;
  double CAMMIN;
  double CAMMAX;
  double LFZO;
  double LCX;
  double LMUX;
  double LEX;
  double LKX;
  double LHX;
  double LVX;
  double LCY;
  double LMUY;
  double LEY;
  double LKY;
  double LHY;
  double LVY;
  double LTR;
  double LRES;
  double LXAL;
  double LYKA;
  double LVYKA;
  double LS;
  double LKYC;
  double LKZC;
  double LVMX;
  double LMX;
  double LMY;
  double LMP;
  double PCX1;
  double PDX1;
  double PDX2;
  double PDX3;
  double PEX1;
  double PEX2;
  double PEX3;
  double PEX4;
  double PKX1;
  double PKX2;
  double PKX3;
  double PHX1;
  double PHX2;
  double PVX1;
  double PVX2;
  double PPX1;
  double PPX2;
  double PPX3;
  double PPX4;
  double RBX1;
  double RBX2;
  double RBX3;
  double RCX1;
  double REX1;
  double REX2;
  double RHX1;
  double QSX1;
  double QSX2;
  double QSX3;
  double QSX4;
  double QSX5;
  double QSX6;
  double QSX7;
  double QSX8;
  double QSX9;
  double QSX10;
  double QSX11;
  double QSX12;
  double QSX13;
  double QSX14;
  double PPMX1;
  double PCY1;
  double PDY1;
  double PDY2;
  double PDY3;
  double PEY1;
  double PEY2;
  double PEY3;
  double PEY4;
  double PEY5;
  double PKY1;
  double PKY2;
  double PKY3;
  double PKY4;
  double PKY5;
  double PKY6;
  double PKY7;
  double PHY1;
  double PHY2;
  double PVY1;
  double PVY2;
  double PVY3;
  double PVY4;
  double PPY1;
  double PPY2;
  double PPY3;
  double PPY4;
  double PPY5;
  double RBY1;
  double RBY2;
  double RBY3;
  double RBY4;
  double RCY1;
  double REY1;
  double REY2;
  double RHY1;
  double RHY2;
  double RVY1;
  double RVY2;
  double RVY3;
  double RVY4;
  double RVY5;
  double RVY6;
  double QSY1;
  double QSY2;
  double QSY3;
  double QSY4;
  double QSY5;
  double QSY6;
  double QSY7;
  double QSY8;
  double QBZ1;
  double QBZ2;
  double QBZ3;
  double QBZ4;
  double QBZ5;
  double QBZ9;
  double QBZ10;
  double QCZ1;
  double QDZ1;
  double QDZ2;
  double QDZ3;
  double QDZ4;
  double QDZ6;
  double QDZ7;
  double QDZ8;
  double QDZ9;
  double QDZ10;
  double QDZ11;
  double QEZ1;
  double QEZ2;
  double QEZ3;
  double QEZ4;
  double QEZ5;
  double QHZ1;
  double QHZ2;
  double QHZ3;
  double QHZ4;
  double PPZ1;
  double PPZ2;
  double SSZ1;
  double SSZ2;
  double SSZ3;
  double SSZ4;
  double PDXP1;
  double PDXP2;
  double PDXP3;
  double PKYP1;
  double PDYP1;
  double PDYP2;
  double PDYP3;
  double PDYP4;
  double PHYP1;
  double PHYP2;
  double PHYP3;
  double PHYP4;
  double PECP1;
  double PECP2;
  double QDTP1;
  double QCRP1;
  double QCRP2;
  double QBRP1;
  double QDRP1;
  char FUNCTION_NAME[39];
  double SWITCH_INTEG;
  double Q_FCY2;
  double Q_CAM1;
  double Q_CAM2;
  double Q_CAM3;
  double Q_FYS1;
  double Q_FYS2;
  double Q_FYS3;
  double ENV_C1;
  double ENV_C2;
  double Q_A1;
  double Q_A2;
  double PHY3;
  double PTX1;
  double PTX2;
  double PTX3;
  double PTY1;
  double PTY2;
  double LSGKP;
  double LSGAL;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  double alpha;
  double kappa;
  double gamma;
  double phit;
  double Fz;
  double p;
  double omega;
  double phi;
  double Vsx;
  double uFz;
  double ukappa;
  double ukappaLow;
  double ualpha;
  double ugamma;
  double uphit;
  double uVcx;
  double nInputs;
  double Fz_lowLimit;
} struct_T;

#endif                                 /*typedef_struct_T*/

#ifndef typedef_cell_0
#define typedef_cell_0

typedef struct {
  struct0_T f2;
  struct_T f3;
  double f4;
  double f5;
  double f6;
  f_struct_T f7;
} cell_0;

#endif                                 /*typedef_cell_0*/

#ifndef struct_emxArray_real_T_1x1
#define struct_emxArray_real_T_1x1

struct emxArray_real_T_1x1
{
  double data[1];
  int size[2];
};

#endif                                 /*struct_emxArray_real_T_1x1*/

#ifndef typedef_emxArray_real_T_1x1
#define typedef_emxArray_real_T_1x1

typedef struct emxArray_real_T_1x1 emxArray_real_T_1x1;

#endif                                 /*typedef_emxArray_real_T_1x1*/

#ifndef typedef_g_struct_T
#define typedef_g_struct_T

typedef struct {
  double Kxk;
  double mux;
  double Kya;
  double muy;
  double t;
  double Mzr;
} g_struct_T;

#endif                                 /*typedef_g_struct_T*/

#ifndef typedef_h_struct_T
#define typedef_h_struct_T

typedef struct {
  double epsilonx;
  double epsilonk;
  double epsilony;
  double epsilonr;
  double epsilonv;
  emxArray_real_T_1x1 reductionSmooth;
  emxArray_real_T_1x1 reductionSharp;
  emxArray_real_T_1x1 reductionLinear;
  emxArray_real_T_1x1 reductionLinear_alpha;
  double zeta2;
  double epsilong;
} h_struct_T;

#endif                                 /*typedef_h_struct_T*/
#endif

/*
 * File trailer for mfeval_wrapper_types.h
 *
 * [EOF]
 */
