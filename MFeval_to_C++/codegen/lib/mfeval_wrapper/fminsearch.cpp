//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fminsearch.cpp
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 09-Jun-2020 11:15:56
//

// Include Files
#include "fminsearch.h"
#include "Solver.h"
#include "mfeval_wrapper.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <math.h>

// Function Definitions

//
// Arguments    : const c_coder_internal_anonymous_func *funfcn
// Return Type  : double
//
double fminsearch(const c_coder_internal_anonymous_func *funfcn)
{
  double x;
  double cfv;
  double xr;
  double fv[2];
  double v[2];
  signed char idx_idx_0;
  signed char idx_idx_1;
  int itercount;
  int fun_evals;
  int lastCol;
  int firstCol;
  boolean_T exitg1;
  int exponent;
  int b_exponent;
  int c_exponent;
  signed char idxb[2];

  //  Cost function to calculate Rl. This measures the error
  //  between the tatget Fz (input of MFeval) vs the Fz calculated
  //  with the current Rl value.
  //  interpolate
  //  Calculate Fz from MF 6.2  with current Rl
  //  Cost (error)
  cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
    (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
    (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6),
    (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy)) -
    funfcn->tunableEnvironment.f7.Fz;
  xr = cfv * cfv;
  fv[0] = xr;
  v[0] = 0.3;
  v[1] = 0.315;

  //  Cost function to calculate Rl. This measures the error
  //  between the tatget Fz (input of MFeval) vs the Fz calculated
  //  with the current Rl value.
  //  interpolate
  //  Calculate Fz from MF 6.2  with current Rl
  //  Cost (error)
  cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
    (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
    (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6), (0.315),
    (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy)) -
    funfcn->tunableEnvironment.f7.Fz;
  cfv *= cfv;
  fv[1] = cfv;
  if ((xr <= cfv) || rtIsNaN(cfv)) {
    idx_idx_0 = 1;
    idx_idx_1 = 2;
  } else {
    idx_idx_0 = 2;
    idx_idx_1 = 1;
  }

  itercount = 1;
  fun_evals = 2;
  lastCol = idx_idx_1 - 1;
  firstCol = idx_idx_0 - 1;
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 200) && (itercount < 200))) {
    int cfv_tmp_tmp;
    int b_cfv_tmp_tmp;
    boolean_T p;
    xr = 0.0;
    cfv_tmp_tmp = idx_idx_0 - 1;
    b_cfv_tmp_tmp = idx_idx_1 - 1;
    cfv = std::abs(fv[cfv_tmp_tmp] - fv[b_cfv_tmp_tmp]);
    if (cfv > 0.0) {
      xr = cfv;
    }

    cfv = std::abs(fv[cfv_tmp_tmp]);
    if ((!rtIsInf(cfv)) && (!rtIsNaN(cfv))) {
      if (cfv <= 2.2250738585072014E-308) {
        cfv = 4.94065645841247E-324;
      } else {
        frexp(cfv, &exponent);
        cfv = std::ldexp(1.0, exponent - 53);
      }
    } else {
      cfv = rtNaN;
    }

    cfv *= 10.0;
    if ((!rtIsInf(cfv)) && (!rtIsNaN(cfv))) {
      if (cfv <= 2.2250738585072014E-308) {
        cfv = 4.94065645841247E-324;
      } else {
        frexp(cfv, &b_exponent);
        cfv = std::ldexp(1.0, b_exponent - 53);
      }
    } else {
      cfv = rtNaN;
    }

    cfv *= 10.0;
    if ((0.0001 > cfv) || rtIsNaN(cfv)) {
      cfv = 0.0001;
    }

    if (xr > cfv) {
      p = false;
    } else {
      xr = 0.0;
      cfv = std::abs(v[cfv_tmp_tmp] - v[b_cfv_tmp_tmp]);
      if (cfv > 0.0) {
        xr = cfv;
      }

      cfv = std::abs(v[cfv_tmp_tmp]);
      if ((!rtIsInf(cfv)) && (!rtIsNaN(cfv))) {
        if (cfv <= 2.2250738585072014E-308) {
          cfv = 4.94065645841247E-324;
        } else {
          frexp(cfv, &c_exponent);
          cfv = std::ldexp(1.0, c_exponent - 53);
        }
      } else {
        cfv = rtNaN;
      }

      cfv *= 10.0;
      if ((1.0E-6 > cfv) || rtIsNaN(cfv)) {
        cfv = 1.0E-6;
      }

      p = (xr <= cfv);
    }

    if (!p) {
      double fxr;
      boolean_T guard1 = false;
      boolean_T guard2 = false;
      xr = 2.0 * v[firstCol] - v[lastCol];

      //  Cost function to calculate Rl. This measures the error
      //  between the tatget Fz (input of MFeval) vs the Fz calculated
      //  with the current Rl value.
      //  interpolate
      //  Calculate Fz from MF 6.2  with current Rl
      //  Cost (error)
      cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
        (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
        (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6), (xr),
        (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy))
        - funfcn->tunableEnvironment.f7.Fz;
      fxr = cfv * cfv;
      fun_evals++;
      guard1 = false;
      guard2 = false;
      if (fxr < fv[cfv_tmp_tmp]) {
        double xe;
        xe = 3.0 * v[firstCol] - 2.0 * v[lastCol];

        //  Cost function to calculate Rl. This measures the error
        //  between the tatget Fz (input of MFeval) vs the Fz calculated
        //  with the current Rl value.
        //  interpolate
        //  Calculate Fz from MF 6.2  with current Rl
        //  Cost (error)
        cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
          (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
          (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6), (xe),
          (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy))
          - funfcn->tunableEnvironment.f7.Fz;
        cfv *= cfv;
        fun_evals++;
        if (cfv < fxr) {
          v[lastCol] = xe;
          fv[b_cfv_tmp_tmp] = cfv;
        } else {
          v[lastCol] = xr;
          fv[b_cfv_tmp_tmp] = fxr;
        }

        guard1 = true;
      } else if (fxr < fv[cfv_tmp_tmp]) {
        v[lastCol] = xr;
        fv[b_cfv_tmp_tmp] = fxr;
        guard1 = true;
      } else if (fxr < fv[b_cfv_tmp_tmp]) {
        x = 1.5 * v[firstCol] - 0.5 * v[lastCol];

        //  Cost function to calculate Rl. This measures the error
        //  between the tatget Fz (input of MFeval) vs the Fz calculated
        //  with the current Rl value.
        //  interpolate
        //  Calculate Fz from MF 6.2  with current Rl
        //  Cost (error)
        cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
          (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
          (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6), (x),
          (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy))
          - funfcn->tunableEnvironment.f7.Fz;
        cfv *= cfv;
        fun_evals++;
        if (cfv <= fxr) {
          v[lastCol] = x;
          fv[b_cfv_tmp_tmp] = cfv;
          guard1 = true;
        } else {
          guard2 = true;
        }
      } else {
        x = 0.5 * v[firstCol] + 0.5 * v[lastCol];

        //  Cost function to calculate Rl. This measures the error
        //  between the tatget Fz (input of MFeval) vs the Fz calculated
        //  with the current Rl value.
        //  interpolate
        //  Calculate Fz from MF 6.2  with current Rl
        //  Cost (error)
        cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
          (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
          (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6), (x),
          (funfcn->tunableEnvironment.f7.Fx), (funfcn->tunableEnvironment.f7.Fy))
          - funfcn->tunableEnvironment.f7.Fz;
        cfv *= cfv;
        fun_evals++;
        if (cfv < fv[b_cfv_tmp_tmp]) {
          v[lastCol] = x;
          fv[b_cfv_tmp_tmp] = cfv;
          guard1 = true;
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        v[b_cfv_tmp_tmp] = v[firstCol] + 0.5 * (v[b_cfv_tmp_tmp] - v[firstCol]);

        //  Cost function to calculate Rl. This measures the error
        //  between the tatget Fz (input of MFeval) vs the Fz calculated
        //  with the current Rl value.
        //  interpolate
        //  Calculate Fz from MF 6.2  with current Rl
        //  Cost (error)
        cfv = mfeval_Solver::calculateFz62((&funfcn->tunableEnvironment.f2),
          (funfcn->tunableEnvironment.f3.gamma), (funfcn->tunableEnvironment.f4),
          (funfcn->tunableEnvironment.f5), (funfcn->tunableEnvironment.f6),
          (v[b_cfv_tmp_tmp]), (funfcn->tunableEnvironment.f7.Fx),
          (funfcn->tunableEnvironment.f7.Fy)) - funfcn->tunableEnvironment.f7.Fz;
        fv[b_cfv_tmp_tmp] = cfv * cfv;
        fun_evals++;
        idxb[0] = idx_idx_0;
        idxb[1] = idx_idx_1;
        if ((fv[cfv_tmp_tmp] <= fv[b_cfv_tmp_tmp]) || rtIsNaN(fv[b_cfv_tmp_tmp]))
        {
          idx_idx_0 = 1;
          idx_idx_1 = 2;
        } else {
          idx_idx_0 = 2;
          idx_idx_1 = 1;
        }

        idx_idx_0 = idxb[idx_idx_0 - 1];
        idx_idx_1 = idxb[idx_idx_1 - 1];
      }

      if (guard1 && (fv[b_cfv_tmp_tmp] < fv[cfv_tmp_tmp])) {
        lastCol = idx_idx_1;
        idx_idx_1 = idx_idx_0;
        idx_idx_0 = static_cast<signed char>(lastCol);
      }

      itercount++;
      lastCol = idx_idx_1 - 1;
      firstCol = idx_idx_0 - 1;
    } else {
      exitg1 = true;
    }
  }

  return v[idx_idx_0 - 1];
}

//
// File trailer for fminsearch.cpp
//
// [EOF]
//
