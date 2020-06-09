/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Solver.h
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

#ifndef SOLVER_H
#define SOLVER_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "mfeval_wrapper_types.h"

/* Function Declarations */
#ifdef __cplusplus

extern "C" {

#endif

  extern void Solver_doForcesAndMoments(const struct0_T *tirParams, const
    struct_T *postProInputs, const h_struct_T *internalParams, const b_struct_T
    modes, f_struct_T *forcesAndmoments, g_struct_T *varinf);
  extern void Solver_parseInputs(const struct0_T *tirParams, const double
    inputs[6], double useMode, struct_T *postProInputs, h_struct_T
    *internalParams, b_struct_T *modes);
  extern double __anon_fcn(const struct0_T *tirParams, const struct_T
    *postProInputs, double omega, double Romega, double dpi, const f_struct_T
    forcesAndmoments, double Rl);

#ifdef __cplusplus

}
#endif
#endif

/*
 * File trailer for Solver.h
 *
 * [EOF]
 */
