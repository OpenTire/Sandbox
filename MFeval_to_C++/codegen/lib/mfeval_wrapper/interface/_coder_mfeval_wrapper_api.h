/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_mfeval_wrapper_api.h
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:15:56
 */

#ifndef _CODER_MFEVAL_WRAPPER_API_H
#define _CODER_MFEVAL_WRAPPER_API_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  char_T FILE_TYPE[5];
  real_T FILE_VERSION;
  char_T FILE_FORMAT[7];
  char_T LENGTH[7];
  char_T FORCE[8];
  char_T ANGLE[9];
  char_T MASS[4];
  char_T TIME[8];
  real_T FITTYP;
  char_T TYRESIDE[6];
  real_T LONGVL;
  real_T VXLOW;
  real_T ROAD_INCREMENT;
  real_T ROAD_DIRECTION;
  char_T PROPERTY_FILE_FORMAT[6];
  real_T USER_SUB_ID;
  real_T N_TIRE_STATES;
  real_T USE_MODE;
  real_T HMAX_LOCAL;
  real_T TIME_SWITCH_INTEG;
  real_T UNLOADED_RADIUS;
  real_T WIDTH;
  real_T ASPECT_RATIO;
  real_T RIM_RADIUS;
  real_T RIM_WIDTH;
  real_T INFLPRES;
  real_T NOMPRES;
  real_T MASS1;
  real_T IXX;
  real_T IYY;
  real_T BELT_MASS;
  real_T BELT_IXX;
  real_T BELT_IYY;
  real_T GRAVITY;
  real_T FNOMIN;
  real_T VERTICAL_STIFFNESS;
  real_T VERTICAL_DAMPING;
  real_T MC_CONTOUR_A;
  real_T MC_CONTOUR_B;
  real_T BREFF;
  real_T DREFF;
  real_T FREFF;
  real_T Q_RE0;
  real_T Q_V1;
  real_T Q_V2;
  real_T Q_FZ2;
  real_T Q_FCX;
  real_T Q_FCY;
  real_T Q_CAM;
  real_T PFZ1;
  real_T BOTTOM_OFFST;
  real_T BOTTOM_STIFF;
  real_T LONGITUDINAL_STIFFNESS;
  real_T LATERAL_STIFFNESS;
  real_T YAW_STIFFNESS;
  real_T FREQ_LONG;
  real_T FREQ_LAT;
  real_T FREQ_YAW;
  real_T FREQ_WINDUP;
  real_T DAMP_LONG;
  real_T DAMP_LAT;
  real_T DAMP_YAW;
  real_T DAMP_WINDUP;
  real_T DAMP_RESIDUAL;
  real_T DAMP_VLOW;
  real_T Q_BVX;
  real_T Q_BVT;
  real_T PCFX1;
  real_T PCFX2;
  real_T PCFX3;
  real_T PCFY1;
  real_T PCFY2;
  real_T PCFY3;
  real_T PCMZ1;
  real_T Q_RA1;
  real_T Q_RA2;
  real_T Q_RB1;
  real_T Q_RB2;
  real_T ELLIPS_SHIFT;
  real_T ELLIPS_LENGTH;
  real_T ELLIPS_HEIGHT;
  real_T ELLIPS_ORDER;
  real_T ELLIPS_MAX_STEP;
  real_T ELLIPS_NWIDTH;
  real_T ELLIPS_NLENGTH;
  real_T PRESMIN;
  real_T PRESMAX;
  real_T FZMIN;
  real_T FZMAX;
  real_T KPUMIN;
  real_T KPUMAX;
  real_T ALPMIN;
  real_T ALPMAX;
  real_T CAMMIN;
  real_T CAMMAX;
  real_T LFZO;
  real_T LCX;
  real_T LMUX;
  real_T LEX;
  real_T LKX;
  real_T LHX;
  real_T LVX;
  real_T LCY;
  real_T LMUY;
  real_T LEY;
  real_T LKY;
  real_T LHY;
  real_T LVY;
  real_T LTR;
  real_T LRES;
  real_T LXAL;
  real_T LYKA;
  real_T LVYKA;
  real_T LS;
  real_T LKYC;
  real_T LKZC;
  real_T LVMX;
  real_T LMX;
  real_T LMY;
  real_T LMP;
  real_T PCX1;
  real_T PDX1;
  real_T PDX2;
  real_T PDX3;
  real_T PEX1;
  real_T PEX2;
  real_T PEX3;
  real_T PEX4;
  real_T PKX1;
  real_T PKX2;
  real_T PKX3;
  real_T PHX1;
  real_T PHX2;
  real_T PVX1;
  real_T PVX2;
  real_T PPX1;
  real_T PPX2;
  real_T PPX3;
  real_T PPX4;
  real_T RBX1;
  real_T RBX2;
  real_T RBX3;
  real_T RCX1;
  real_T REX1;
  real_T REX2;
  real_T RHX1;
  real_T QSX1;
  real_T QSX2;
  real_T QSX3;
  real_T QSX4;
  real_T QSX5;
  real_T QSX6;
  real_T QSX7;
  real_T QSX8;
  real_T QSX9;
  real_T QSX10;
  real_T QSX11;
  real_T QSX12;
  real_T QSX13;
  real_T QSX14;
  real_T PPMX1;
  real_T PCY1;
  real_T PDY1;
  real_T PDY2;
  real_T PDY3;
  real_T PEY1;
  real_T PEY2;
  real_T PEY3;
  real_T PEY4;
  real_T PEY5;
  real_T PKY1;
  real_T PKY2;
  real_T PKY3;
  real_T PKY4;
  real_T PKY5;
  real_T PKY6;
  real_T PKY7;
  real_T PHY1;
  real_T PHY2;
  real_T PVY1;
  real_T PVY2;
  real_T PVY3;
  real_T PVY4;
  real_T PPY1;
  real_T PPY2;
  real_T PPY3;
  real_T PPY4;
  real_T PPY5;
  real_T RBY1;
  real_T RBY2;
  real_T RBY3;
  real_T RBY4;
  real_T RCY1;
  real_T REY1;
  real_T REY2;
  real_T RHY1;
  real_T RHY2;
  real_T RVY1;
  real_T RVY2;
  real_T RVY3;
  real_T RVY4;
  real_T RVY5;
  real_T RVY6;
  real_T QSY1;
  real_T QSY2;
  real_T QSY3;
  real_T QSY4;
  real_T QSY5;
  real_T QSY6;
  real_T QSY7;
  real_T QSY8;
  real_T QBZ1;
  real_T QBZ2;
  real_T QBZ3;
  real_T QBZ4;
  real_T QBZ5;
  real_T QBZ9;
  real_T QBZ10;
  real_T QCZ1;
  real_T QDZ1;
  real_T QDZ2;
  real_T QDZ3;
  real_T QDZ4;
  real_T QDZ6;
  real_T QDZ7;
  real_T QDZ8;
  real_T QDZ9;
  real_T QDZ10;
  real_T QDZ11;
  real_T QEZ1;
  real_T QEZ2;
  real_T QEZ3;
  real_T QEZ4;
  real_T QEZ5;
  real_T QHZ1;
  real_T QHZ2;
  real_T QHZ3;
  real_T QHZ4;
  real_T PPZ1;
  real_T PPZ2;
  real_T SSZ1;
  real_T SSZ2;
  real_T SSZ3;
  real_T SSZ4;
  real_T PDXP1;
  real_T PDXP2;
  real_T PDXP3;
  real_T PKYP1;
  real_T PDYP1;
  real_T PDYP2;
  real_T PDYP3;
  real_T PDYP4;
  real_T PHYP1;
  real_T PHYP2;
  real_T PHYP3;
  real_T PHYP4;
  real_T PECP1;
  real_T PECP2;
  real_T QDTP1;
  real_T QCRP1;
  real_T QCRP2;
  real_T QBRP1;
  real_T QDRP1;
  char_T FUNCTION_NAME[39];
  real_T SWITCH_INTEG;
  real_T Q_FCY2;
  real_T Q_CAM1;
  real_T Q_CAM2;
  real_T Q_CAM3;
  real_T Q_FYS1;
  real_T Q_FYS2;
  real_T Q_FYS3;
  real_T ENV_C1;
  real_T ENV_C2;
  real_T Q_A1;
  real_T Q_A2;
  real_T PHY3;
  real_T PTX1;
  real_T PTX2;
  real_T PTX3;
  real_T PTY1;
  real_T PTY2;
  real_T LSGKP;
  real_T LSGAL;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void mfeval_wrapper(struct0_T *parameterSource, real_T Fz, real_T kappa,
  real_T alpha, real_T b_gamma, real_T phit, real_T Vx, real_T useMode, real_T
  outMF[30]);
extern void mfeval_wrapper_api(const mxArray * const prhs[8], int32_T nlhs,
  const mxArray *plhs[1]);
extern void mfeval_wrapper_atexit(void);
extern void mfeval_wrapper_initialize(void);
extern void mfeval_wrapper_terminate(void);
extern void mfeval_wrapper_xil_shutdown(void);
extern void mfeval_wrapper_xil_terminate(void);

#endif

/*
 * File trailer for _coder_mfeval_wrapper_api.h
 *
 * [EOF]
 */
