//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: main.cpp
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 09-Jun-2020 11:15:56
//

//***********************************************************************
// This automatically generated example C++ main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************

// Include Files
#include "main.h"
#include "mfeval_wrapper.h"
#include "mfeval_wrapper_terminate.h"
#include "rt_nonfinite.h"

// Function Declarations
static void argInit_1x39_char_T(char result[39]);
static void argInit_1x4_char_T(char result[4]);
static void argInit_1x5_char_T(char result[5]);
static void argInit_1x6_char_T(char result[6]);
static void argInit_1x7_char_T(char result[7]);
static void argInit_1x8_char_T(char result[8]);
static void argInit_1x9_char_T(char result[9]);
static char argInit_char_T();
static double argInit_real_T();
static void argInit_struct0_T(struct0_T *result);
static void main_mfeval_wrapper();

// Function Definitions

//
// Arguments    : char result[39]
// Return Type  : void
//
static void argInit_1x39_char_T(char result[39])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 39; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : char result[4]
// Return Type  : void
//
static void argInit_1x4_char_T(char result[4])
{
  char result_tmp;

  // Loop over the array to initialize each element.
  // Set the value of the array element.
  // Change this value to the value that the application requires.
  result_tmp = argInit_char_T();
  result[0] = result_tmp;

  // Set the value of the array element.
  // Change this value to the value that the application requires.
  result[1] = result_tmp;

  // Set the value of the array element.
  // Change this value to the value that the application requires.
  result[2] = result_tmp;

  // Set the value of the array element.
  // Change this value to the value that the application requires.
  result[3] = result_tmp;
}

//
// Arguments    : char result[5]
// Return Type  : void
//
static void argInit_1x5_char_T(char result[5])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 5; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : char result[6]
// Return Type  : void
//
static void argInit_1x6_char_T(char result[6])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 6; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : char result[7]
// Return Type  : void
//
static void argInit_1x7_char_T(char result[7])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 7; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : char result[8]
// Return Type  : void
//
static void argInit_1x8_char_T(char result[8])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 8; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : char result[9]
// Return Type  : void
//
static void argInit_1x9_char_T(char result[9])
{
  // Loop over the array to initialize each element.
  for (int idx1 = 0; idx1 < 9; idx1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx1] = argInit_char_T();
  }
}

//
// Arguments    : void
// Return Type  : char
//
static char argInit_char_T()
{
  return '?';
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : struct0_T *result
// Return Type  : void
//
static void argInit_struct0_T(struct0_T *result)
{
  char result_tmp[7];
  char b_result_tmp[8];
  double c_result_tmp;
  char d_result_tmp[6];
  int i;

  // Set the value of each structure field.
  // Change this value to the value that the application requires.
  argInit_1x7_char_T(result_tmp);
  argInit_1x8_char_T(b_result_tmp);
  c_result_tmp = argInit_real_T();
  result->FITTYP = c_result_tmp;
  result->LONGVL = c_result_tmp;
  result->VXLOW = c_result_tmp;
  result->ROAD_INCREMENT = c_result_tmp;
  result->ROAD_DIRECTION = c_result_tmp;
  argInit_1x6_char_T(d_result_tmp);
  result->USER_SUB_ID = c_result_tmp;
  result->N_TIRE_STATES = c_result_tmp;
  result->USE_MODE = c_result_tmp;
  result->HMAX_LOCAL = c_result_tmp;
  result->TIME_SWITCH_INTEG = c_result_tmp;
  result->UNLOADED_RADIUS = c_result_tmp;
  result->WIDTH = c_result_tmp;
  result->ASPECT_RATIO = c_result_tmp;
  result->RIM_RADIUS = c_result_tmp;
  result->RIM_WIDTH = c_result_tmp;
  result->INFLPRES = c_result_tmp;
  result->NOMPRES = c_result_tmp;
  result->MASS1 = c_result_tmp;
  result->IXX = c_result_tmp;
  result->IYY = c_result_tmp;
  result->BELT_MASS = c_result_tmp;
  result->BELT_IXX = c_result_tmp;
  result->BELT_IYY = c_result_tmp;
  result->GRAVITY = c_result_tmp;
  result->FNOMIN = c_result_tmp;
  result->VERTICAL_STIFFNESS = c_result_tmp;
  result->VERTICAL_DAMPING = c_result_tmp;
  result->MC_CONTOUR_A = c_result_tmp;
  result->MC_CONTOUR_B = c_result_tmp;
  result->BREFF = c_result_tmp;
  result->DREFF = c_result_tmp;
  result->FREFF = c_result_tmp;
  result->Q_RE0 = c_result_tmp;
  result->Q_V1 = c_result_tmp;
  result->Q_V2 = c_result_tmp;
  result->Q_FZ2 = c_result_tmp;
  result->Q_FCX = c_result_tmp;
  result->Q_FCY = c_result_tmp;
  result->Q_CAM = c_result_tmp;
  result->PFZ1 = c_result_tmp;
  result->BOTTOM_OFFST = c_result_tmp;
  result->BOTTOM_STIFF = c_result_tmp;
  result->LONGITUDINAL_STIFFNESS = c_result_tmp;
  result->LATERAL_STIFFNESS = c_result_tmp;
  result->YAW_STIFFNESS = c_result_tmp;
  result->FREQ_LONG = c_result_tmp;
  result->FREQ_LAT = c_result_tmp;
  result->FREQ_YAW = c_result_tmp;
  result->FREQ_WINDUP = c_result_tmp;
  result->DAMP_LONG = c_result_tmp;
  result->DAMP_LAT = c_result_tmp;
  result->DAMP_YAW = c_result_tmp;
  result->DAMP_WINDUP = c_result_tmp;
  result->DAMP_RESIDUAL = c_result_tmp;
  result->DAMP_VLOW = c_result_tmp;
  result->Q_BVX = c_result_tmp;
  result->Q_BVT = c_result_tmp;
  result->PCFX1 = c_result_tmp;
  result->PCFX2 = c_result_tmp;
  result->PCFX3 = c_result_tmp;
  result->PCFY1 = c_result_tmp;
  result->PCFY2 = c_result_tmp;
  result->PCFY3 = c_result_tmp;
  result->PCMZ1 = c_result_tmp;
  result->Q_RA1 = c_result_tmp;
  result->Q_RA2 = c_result_tmp;
  result->Q_RB1 = c_result_tmp;
  result->Q_RB2 = c_result_tmp;
  result->ELLIPS_SHIFT = c_result_tmp;
  result->ELLIPS_LENGTH = c_result_tmp;
  result->ELLIPS_HEIGHT = c_result_tmp;
  result->ELLIPS_ORDER = c_result_tmp;
  result->ELLIPS_MAX_STEP = c_result_tmp;
  result->ELLIPS_NWIDTH = c_result_tmp;
  result->ELLIPS_NLENGTH = c_result_tmp;
  result->PRESMIN = c_result_tmp;
  result->PRESMAX = c_result_tmp;
  result->FZMIN = c_result_tmp;
  result->FZMAX = c_result_tmp;
  result->KPUMIN = c_result_tmp;
  result->KPUMAX = c_result_tmp;
  result->ALPMIN = c_result_tmp;
  result->ALPMAX = c_result_tmp;
  result->CAMMIN = c_result_tmp;
  result->CAMMAX = c_result_tmp;
  result->LFZO = c_result_tmp;
  result->LCX = c_result_tmp;
  result->LMUX = c_result_tmp;
  result->LEX = c_result_tmp;
  result->LKX = c_result_tmp;
  result->LHX = c_result_tmp;
  result->LVX = c_result_tmp;
  result->LCY = c_result_tmp;
  result->LMUY = c_result_tmp;
  result->LEY = c_result_tmp;
  result->LKY = c_result_tmp;
  result->LHY = c_result_tmp;
  result->LVY = c_result_tmp;
  result->LTR = c_result_tmp;
  result->LRES = c_result_tmp;
  result->LXAL = c_result_tmp;
  result->LYKA = c_result_tmp;
  result->LVYKA = c_result_tmp;
  result->LS = c_result_tmp;
  result->LKYC = c_result_tmp;
  result->LKZC = c_result_tmp;
  result->LVMX = c_result_tmp;
  result->LMX = c_result_tmp;
  result->LMY = c_result_tmp;
  result->LMP = c_result_tmp;
  result->PCX1 = c_result_tmp;
  result->PDX1 = c_result_tmp;
  result->PDX2 = c_result_tmp;
  result->PDX3 = c_result_tmp;
  result->PEX1 = c_result_tmp;
  result->PEX2 = c_result_tmp;
  result->PEX3 = c_result_tmp;
  result->PEX4 = c_result_tmp;
  result->PKX1 = c_result_tmp;
  result->PKX2 = c_result_tmp;
  result->PKX3 = c_result_tmp;
  result->PHX1 = c_result_tmp;
  result->PHX2 = c_result_tmp;
  result->PVX1 = c_result_tmp;
  result->PVX2 = c_result_tmp;
  result->PPX1 = c_result_tmp;
  result->PPX2 = c_result_tmp;
  result->PPX3 = c_result_tmp;
  result->PPX4 = c_result_tmp;
  result->RBX1 = c_result_tmp;
  result->RBX2 = c_result_tmp;
  result->RBX3 = c_result_tmp;
  result->RCX1 = c_result_tmp;
  result->REX1 = c_result_tmp;
  result->REX2 = c_result_tmp;
  result->RHX1 = c_result_tmp;
  result->QSX1 = c_result_tmp;
  result->QSX2 = c_result_tmp;
  result->QSX3 = c_result_tmp;
  result->QSX4 = c_result_tmp;
  result->QSX5 = c_result_tmp;
  result->QSX6 = c_result_tmp;
  result->QSX7 = c_result_tmp;
  result->QSX8 = c_result_tmp;
  result->QSX9 = c_result_tmp;
  result->QSX10 = c_result_tmp;
  result->QSX11 = c_result_tmp;
  result->QSX12 = c_result_tmp;
  result->QSX13 = c_result_tmp;
  result->QSX14 = c_result_tmp;
  result->PPMX1 = c_result_tmp;
  result->PCY1 = c_result_tmp;
  result->PDY1 = c_result_tmp;
  result->PDY2 = c_result_tmp;
  result->PDY3 = c_result_tmp;
  result->PEY1 = c_result_tmp;
  result->PEY2 = c_result_tmp;
  result->PEY3 = c_result_tmp;
  result->PEY4 = c_result_tmp;
  result->PEY5 = c_result_tmp;
  result->PKY1 = c_result_tmp;
  result->PKY2 = c_result_tmp;
  result->PKY3 = c_result_tmp;
  result->PKY4 = c_result_tmp;
  result->PKY5 = c_result_tmp;
  result->PKY6 = c_result_tmp;
  result->PKY7 = c_result_tmp;
  result->PHY1 = c_result_tmp;
  result->PHY2 = c_result_tmp;
  result->PVY1 = c_result_tmp;
  result->PVY2 = c_result_tmp;
  result->PVY3 = c_result_tmp;
  result->PVY4 = c_result_tmp;
  result->PPY1 = c_result_tmp;
  result->PPY2 = c_result_tmp;
  result->PPY3 = c_result_tmp;
  result->PPY4 = c_result_tmp;
  result->PPY5 = c_result_tmp;
  result->RBY1 = c_result_tmp;
  result->RBY2 = c_result_tmp;
  result->RBY3 = c_result_tmp;
  result->RBY4 = c_result_tmp;
  result->RCY1 = c_result_tmp;
  result->REY1 = c_result_tmp;
  result->REY2 = c_result_tmp;
  result->RHY1 = c_result_tmp;
  result->RHY2 = c_result_tmp;
  result->RVY1 = c_result_tmp;
  result->RVY2 = c_result_tmp;
  result->RVY3 = c_result_tmp;
  result->RVY4 = c_result_tmp;
  result->RVY5 = c_result_tmp;
  result->RVY6 = c_result_tmp;
  result->QSY1 = c_result_tmp;
  result->QSY2 = c_result_tmp;
  result->QSY3 = c_result_tmp;
  result->QSY4 = c_result_tmp;
  result->QSY5 = c_result_tmp;
  result->QSY6 = c_result_tmp;
  result->QSY7 = c_result_tmp;
  result->QSY8 = c_result_tmp;
  result->QBZ1 = c_result_tmp;
  result->QBZ2 = c_result_tmp;
  result->QBZ3 = c_result_tmp;
  result->QBZ4 = c_result_tmp;
  result->QBZ5 = c_result_tmp;
  result->QBZ9 = c_result_tmp;
  result->QBZ10 = c_result_tmp;
  result->QCZ1 = c_result_tmp;
  result->QDZ1 = c_result_tmp;
  result->QDZ2 = c_result_tmp;
  result->QDZ3 = c_result_tmp;
  result->QDZ4 = c_result_tmp;
  result->QDZ6 = c_result_tmp;
  result->QDZ7 = c_result_tmp;
  result->QDZ8 = c_result_tmp;
  result->QDZ9 = c_result_tmp;
  result->QDZ10 = c_result_tmp;
  result->QDZ11 = c_result_tmp;
  result->QEZ1 = c_result_tmp;
  result->QEZ2 = c_result_tmp;
  result->QEZ3 = c_result_tmp;
  result->QEZ4 = c_result_tmp;
  result->QEZ5 = c_result_tmp;
  result->QHZ1 = c_result_tmp;
  result->QHZ2 = c_result_tmp;
  result->QHZ3 = c_result_tmp;
  result->QHZ4 = c_result_tmp;
  result->PPZ1 = c_result_tmp;
  result->PPZ2 = c_result_tmp;
  result->SSZ1 = c_result_tmp;
  result->SSZ2 = c_result_tmp;
  result->SSZ3 = c_result_tmp;
  result->SSZ4 = c_result_tmp;
  result->PDXP1 = c_result_tmp;
  result->PDXP2 = c_result_tmp;
  result->PDXP3 = c_result_tmp;
  result->PKYP1 = c_result_tmp;
  result->PDYP1 = c_result_tmp;
  result->PDYP2 = c_result_tmp;
  result->PDYP3 = c_result_tmp;
  result->PDYP4 = c_result_tmp;
  result->PHYP1 = c_result_tmp;
  result->PHYP2 = c_result_tmp;
  result->PHYP3 = c_result_tmp;
  result->PHYP4 = c_result_tmp;
  result->PECP1 = c_result_tmp;
  result->PECP2 = c_result_tmp;
  result->QDTP1 = c_result_tmp;
  result->QCRP1 = c_result_tmp;
  result->QCRP2 = c_result_tmp;
  result->QBRP1 = c_result_tmp;
  result->QDRP1 = c_result_tmp;
  result->SWITCH_INTEG = c_result_tmp;
  result->Q_FCY2 = c_result_tmp;
  result->Q_CAM1 = c_result_tmp;
  result->Q_CAM2 = c_result_tmp;
  result->Q_CAM3 = c_result_tmp;
  result->Q_FYS1 = c_result_tmp;
  result->Q_FYS2 = c_result_tmp;
  result->Q_FYS3 = c_result_tmp;
  result->ENV_C1 = c_result_tmp;
  result->ENV_C2 = c_result_tmp;
  result->Q_A1 = c_result_tmp;
  result->Q_A2 = c_result_tmp;
  result->PHY3 = c_result_tmp;
  result->PTX1 = c_result_tmp;
  result->PTX2 = c_result_tmp;
  result->PTX3 = c_result_tmp;
  result->PTY1 = c_result_tmp;
  result->PTY2 = c_result_tmp;
  result->LSGKP = c_result_tmp;
  result->LSGAL = c_result_tmp;
  argInit_1x5_char_T(result->FILE_TYPE);
  result->FILE_VERSION = c_result_tmp;
  argInit_1x9_char_T(result->ANGLE);
  argInit_1x4_char_T(result->MASS);
  argInit_1x39_char_T(result->FUNCTION_NAME);
  for (i = 0; i < 7; i++) {
    result->FILE_FORMAT[i] = result_tmp[i];
    result->LENGTH[i] = result_tmp[i];
  }

  for (i = 0; i < 8; i++) {
    result->FORCE[i] = b_result_tmp[i];
    result->TIME[i] = b_result_tmp[i];
  }

  for (i = 0; i < 6; i++) {
    result->TYRESIDE[i] = d_result_tmp[i];
    result->PROPERTY_FILE_FORMAT[i] = d_result_tmp[i];
  }
}

//
// Arguments    : void
// Return Type  : void
//
static void main_mfeval_wrapper()
{
  double Fz_tmp;
  struct0_T r;
  double outMF[30];

  // Initialize function 'mfeval_wrapper' input arguments.
  // Initialize function input argument 'parameterSource'.
  Fz_tmp = argInit_real_T();

  // Call the entry-point 'mfeval_wrapper'.
  argInit_struct0_T(&r);
  mfeval_wrapper(&r, Fz_tmp, Fz_tmp, Fz_tmp, Fz_tmp, Fz_tmp, Fz_tmp, Fz_tmp,
                 outMF);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. 
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_mfeval_wrapper();

  // Terminate the application.
  // You do not need to do this more than one time.
  mfeval_wrapper_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
