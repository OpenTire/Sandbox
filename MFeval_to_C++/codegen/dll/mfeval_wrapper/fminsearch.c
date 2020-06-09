/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: fminsearch.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

/* Include Files */
#include "fminsearch.h"
#include "Solver.h"
#include "mfeval_wrapper.h"
#include "mfeval_wrapper_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */

/*
 * Arguments    : const cell_0 *funfcn_tunableEnvironment
 * Return Type  : double
 */
double fminsearch(const cell_0 *funfcn_tunableEnvironment)
{
  double x;
  double a;
  double a_tmp;
  double rho_zfr;
  double fxr;
  double fxe;
  double rho_zg;
  double dv[1];
  int lastCol;
  double fv[2];
  double v[2];
  signed char idx_idx_0;
  signed char idx_idx_1;
  int itercount;
  int fun_evals;
  int firstCol;
  boolean_T exitg1;
  int cfv_tmp_tmp;
  int b_cfv_tmp_tmp;
  int exponent;
  int b_exponent;
  boolean_T p;
  boolean_T guard1 = false;
  int c_exponent;
  boolean_T guard2 = false;
  signed char idxb[2];

  /*  Cost function to calculate Rl. This measures the error */
  /*  between the tatget Fz (input of MFeval) vs the Fz calculated */
  /*  with the current Rl value. */
  /*  interpolate */
  /*  Calculate Fz from MF 6.2  with current Rl */
  /*  Cost (error) */
  a = funfcn_tunableEnvironment->f2.VERTICAL_STIFFNESS *
    funfcn_tunableEnvironment->f2.UNLOADED_RADIUS /
    funfcn_tunableEnvironment->f2.FNOMIN;
  a_tmp = 0.3 / funfcn_tunableEnvironment->f5;
  rho_zfr = funfcn_tunableEnvironment->f5 - 0.3;
  if (!(rho_zfr > 0.0)) {
    rho_zfr = 0.0;
  }

  fxr = (funfcn_tunableEnvironment->f2.Q_CAM1 * 0.3 +
         funfcn_tunableEnvironment->f2.Q_CAM2 * 0.09) *
    funfcn_tunableEnvironment->f3.gamma;
  fxe = (funfcn_tunableEnvironment->f2.Q_CAM1 * funfcn_tunableEnvironment->f5 +
         funfcn_tunableEnvironment->f2.Q_CAM2 * (funfcn_tunableEnvironment->f5 *
          funfcn_tunableEnvironment->f5)) * funfcn_tunableEnvironment->f3.gamma;
  rho_zg = fxr * fxr * ((1.075 - 0.5 *
    funfcn_tunableEnvironment->f2.ASPECT_RATIO) *
                        funfcn_tunableEnvironment->f2.WIDTH / 8.0) * fabs(tan
    (funfcn_tunableEnvironment->f3.gamma)) / (fxe * fxe) -
    funfcn_tunableEnvironment->f2.Q_CAM3 * rho_zfr * fabs
    (funfcn_tunableEnvironment->f3.gamma);
  dv[0] = rho_zg;
  lastCol = 0;
  if (rtIsNaN(rho_zg)) {
    lastCol = 1;
  }

  if (0 <= lastCol - 1) {
    dv[0] = 0.0;
  }

  fxr = funfcn_tunableEnvironment->f2.Q_FCX * funfcn_tunableEnvironment->f7.Fx /
    funfcn_tunableEnvironment->f2.FNOMIN;
  rho_zg = rho_zfr + dv[0];
  if (!(rho_zg > 0.0)) {
    rho_zg = 0.0;
  }

  rho_zg /= funfcn_tunableEnvironment->f2.UNLOADED_RADIUS;
  fxe = rt_powd_snf(rho_zg, funfcn_tunableEnvironment->f2.Q_FCY2) *
    (funfcn_tunableEnvironment->f2.Q_FCY * (funfcn_tunableEnvironment->f7.Fy -
      ((funfcn_tunableEnvironment->f2.Q_FYS1 +
        funfcn_tunableEnvironment->f2.Q_FYS2 * a_tmp) +
       funfcn_tunableEnvironment->f2.Q_FYS3 * (a_tmp * a_tmp)) *
      funfcn_tunableEnvironment->f3.gamma) /
     funfcn_tunableEnvironment->f2.FNOMIN);
  a = (((funfcn_tunableEnvironment->f2.Q_V2 *
         (funfcn_tunableEnvironment->f2.UNLOADED_RADIUS /
          funfcn_tunableEnvironment->f2.LONGVL) * fabs
         (funfcn_tunableEnvironment->f4) + 1.0) - fxr * fxr) - fxe * fxe) *
    (funfcn_tunableEnvironment->f2.PFZ1 * funfcn_tunableEnvironment->f6 + 1.0) *
    (sqrt(a * a - 4.0 * funfcn_tunableEnvironment->f2.Q_FZ2) * rho_zg +
     funfcn_tunableEnvironment->f2.Q_FZ2 * (rho_zg * rho_zg)) *
    funfcn_tunableEnvironment->f2.FNOMIN - funfcn_tunableEnvironment->f7.Fz;
  rho_zg = a * a;
  fv[0] = rho_zg;
  v[0] = 0.3;
  v[1] = 0.315;
  fv[1] = __anon_fcn(&funfcn_tunableEnvironment->f2,
                     &funfcn_tunableEnvironment->f3,
                     funfcn_tunableEnvironment->f4,
                     funfcn_tunableEnvironment->f5,
                     funfcn_tunableEnvironment->f6,
                     funfcn_tunableEnvironment->f7, 0.315);
  if ((rho_zg <= fv[1]) || rtIsNaN(fv[1])) {
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
    rho_zfr = 0.0;
    cfv_tmp_tmp = idx_idx_0 - 1;
    b_cfv_tmp_tmp = idx_idx_1 - 1;
    rho_zg = fabs(fv[cfv_tmp_tmp] - fv[b_cfv_tmp_tmp]);
    if (rho_zg > 0.0) {
      rho_zfr = rho_zg;
    }

    rho_zg = fabs(fv[cfv_tmp_tmp]);
    if ((!rtIsInf(rho_zg)) && (!rtIsNaN(rho_zg))) {
      if (rho_zg <= 2.2250738585072014E-308) {
        rho_zg = 4.94065645841247E-324;
      } else {
        frexp(rho_zg, &exponent);
        rho_zg = ldexp(1.0, exponent - 53);
      }
    } else {
      rho_zg = rtNaN;
    }

    rho_zg *= 10.0;
    if ((!rtIsInf(rho_zg)) && (!rtIsNaN(rho_zg))) {
      if (rho_zg <= 2.2250738585072014E-308) {
        rho_zg = 4.94065645841247E-324;
      } else {
        frexp(rho_zg, &b_exponent);
        rho_zg = ldexp(1.0, b_exponent - 53);
      }
    } else {
      rho_zg = rtNaN;
    }

    rho_zg *= 10.0;
    if ((0.0001 > rho_zg) || rtIsNaN(rho_zg)) {
      rho_zg = 0.0001;
    }

    if (rho_zfr > rho_zg) {
      p = false;
    } else {
      rho_zfr = 0.0;
      rho_zg = fabs(v[cfv_tmp_tmp] - v[b_cfv_tmp_tmp]);
      if (rho_zg > 0.0) {
        rho_zfr = rho_zg;
      }

      rho_zg = fabs(v[cfv_tmp_tmp]);
      if ((!rtIsInf(rho_zg)) && (!rtIsNaN(rho_zg))) {
        if (rho_zg <= 2.2250738585072014E-308) {
          rho_zg = 4.94065645841247E-324;
        } else {
          frexp(rho_zg, &c_exponent);
          rho_zg = ldexp(1.0, c_exponent - 53);
        }
      } else {
        rho_zg = rtNaN;
      }

      rho_zg *= 10.0;
      if ((1.0E-6 > rho_zg) || rtIsNaN(rho_zg)) {
        rho_zg = 1.0E-6;
      }

      p = (rho_zfr <= rho_zg);
    }

    if (!p) {
      rho_zg = 2.0 * v[firstCol] - v[lastCol];
      fxr = __anon_fcn(&funfcn_tunableEnvironment->f2,
                       &funfcn_tunableEnvironment->f3,
                       funfcn_tunableEnvironment->f4,
                       funfcn_tunableEnvironment->f5,
                       funfcn_tunableEnvironment->f6,
                       funfcn_tunableEnvironment->f7, rho_zg);
      fun_evals++;
      guard1 = false;
      guard2 = false;
      if (fxr < fv[cfv_tmp_tmp]) {
        rho_zfr = 3.0 * v[firstCol] - 2.0 * v[lastCol];
        fxe = __anon_fcn(&funfcn_tunableEnvironment->f2,
                         &funfcn_tunableEnvironment->f3,
                         funfcn_tunableEnvironment->f4,
                         funfcn_tunableEnvironment->f5,
                         funfcn_tunableEnvironment->f6,
                         funfcn_tunableEnvironment->f7, rho_zfr);
        fun_evals++;
        if (fxe < fxr) {
          v[lastCol] = rho_zfr;
          fv[b_cfv_tmp_tmp] = fxe;
        } else {
          v[lastCol] = rho_zg;
          fv[b_cfv_tmp_tmp] = fxr;
        }

        guard1 = true;
      } else if (fxr < fv[cfv_tmp_tmp]) {
        v[lastCol] = rho_zg;
        fv[b_cfv_tmp_tmp] = fxr;
        guard1 = true;
      } else if (fxr < fv[b_cfv_tmp_tmp]) {
        x = 1.5 * v[firstCol] - 0.5 * v[lastCol];
        rho_zg = __anon_fcn(&funfcn_tunableEnvironment->f2,
                            &funfcn_tunableEnvironment->f3,
                            funfcn_tunableEnvironment->f4,
                            funfcn_tunableEnvironment->f5,
                            funfcn_tunableEnvironment->f6,
                            funfcn_tunableEnvironment->f7, x);
        fun_evals++;
        if (rho_zg <= fxr) {
          v[lastCol] = x;
          fv[b_cfv_tmp_tmp] = rho_zg;
          guard1 = true;
        } else {
          guard2 = true;
        }
      } else {
        x = 0.5 * v[firstCol] + 0.5 * v[lastCol];
        rho_zg = __anon_fcn(&funfcn_tunableEnvironment->f2,
                            &funfcn_tunableEnvironment->f3,
                            funfcn_tunableEnvironment->f4,
                            funfcn_tunableEnvironment->f5,
                            funfcn_tunableEnvironment->f6,
                            funfcn_tunableEnvironment->f7, x);
        fun_evals++;
        if (rho_zg < fv[b_cfv_tmp_tmp]) {
          v[lastCol] = x;
          fv[b_cfv_tmp_tmp] = rho_zg;
          guard1 = true;
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        v[b_cfv_tmp_tmp] = v[firstCol] + 0.5 * (v[b_cfv_tmp_tmp] - v[firstCol]);
        fv[b_cfv_tmp_tmp] = __anon_fcn(&funfcn_tunableEnvironment->f2,
          &funfcn_tunableEnvironment->f3, funfcn_tunableEnvironment->f4,
          funfcn_tunableEnvironment->f5, funfcn_tunableEnvironment->f6,
          funfcn_tunableEnvironment->f7, v[b_cfv_tmp_tmp]);
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
        idx_idx_0 = (signed char)lastCol;
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

/*
 * File trailer for fminsearch.c
 *
 * [EOF]
 */
