/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mfeval_wrapper_initialize.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

/* Include Files */
#include "mfeval_wrapper_initialize.h"
#include "mfeval_wrapper.h"
#include "mfeval_wrapper_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void mfeval_wrapper_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_mfeval_wrapper = true;
}

/*
 * File trailer for mfeval_wrapper_initialize.c
 *
 * [EOF]
 */
