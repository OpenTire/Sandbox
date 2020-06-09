/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Solver.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

/* Include Files */
#include "Solver.h"
#include "mfeval_wrapper.h"
#include "mfeval_wrapper_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_c_struct_T
#define typedef_c_struct_T

typedef struct {
  double dfz;
  double dpi;
} c_struct_T;

#endif                                 /*typedef_c_struct_T*/

#ifndef typedef_d_struct_T
#define typedef_d_struct_T

typedef struct {
  double Fz0_prime;
  double alpha_prime;
  double LMUX_prime;
  double LMUY_prime;
} d_struct_T;

#endif                                 /*typedef_d_struct_T*/

#ifndef typedef_e_struct_T
#define typedef_e_struct_T

typedef struct {
  double alpha_star;
  double gamma_star;
  double LMUX_star;
  double LMUY_star;
} e_struct_T;

#endif                                 /*typedef_e_struct_T*/

#ifndef typedef_i_struct_T
#define typedef_i_struct_T

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
  double zeta0;
  double zeta1;
  double zeta3;
  double zeta4;
  double zeta5;
  double zeta6;
  double zeta7;
  double zeta8;
} i_struct_T;

#endif                                 /*typedef_i_struct_T*/

/* Function Declarations */
static void Solver_calculateFx0(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fx0, double *mux, double *Kxk);
static void Solver_calculateFy0(const struct0_T *tirParams, const struct_T
  *postProInputs, const i_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fy0, double *muy, double *Kya, double *Kyg0, double *SHy, double *SVy,
  double *By, double *Cy);
static void Solver_calculateMz(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double alphar, double alphat, double Kxk, double Kya_prime, double Fy, double
  Fx, double Dr, double Cr, double Br, double Dt, double Ct, double Bt, double
  Et, double *Mz, double *t, double *Mzr);
static void Solver_calculateMz0(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double Kya, double SHy, double SVy, double By, double Cy, double *Mz0, double *
  alphar, double *alphat, double *Dr, double *Cr, double *Br, double *Dt, double
  *Ct, double *Bt, double *Et, double *Kya_prime);
static void Solver_calculateRe(double tirParams_LONGVL, double
  tirParams_UNLOADED_RADIUS, double tirParams_FNOMIN, double
  tirParams_VERTICAL_STIFFNESS, double tirParams_BREFF, double tirParams_DREFF,
  double tirParams_FREFF, double tirParams_Q_RE0, double tirParams_Q_V1, double
  tirParams_PFZ1, double postProInputs_omega, double postProInputs_uFz, double
  postProInputs_ukappa, double postProInputs_uVcx, double dpi, double *Re,
  double *Romega, double *omega);
static void b_Solver_calculateFy0(const struct0_T *tirParams, const struct_T
  *postProInputs, h_struct_T *internalParams, const b_struct_T modes, const
  e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fy0, double *muy, double *Kya, double *Kyg0, double *SHy, double *SVy,
  double *By, double *Cy);

/* Function Definitions */

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                const h_struct_T *internalParams
 *                const b_struct_T modes
 *                const e_struct_T starVar
 *                const d_struct_T primeVar
 *                const c_struct_T incrVar
 *                double *Fx0
 *                double *mux
 *                double *Kxk
 * Return Type  : void
 */
static void Solver_calculateFx0(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fx0, double *mux, double *Kxk)
{
  double zeta1;
  double Cx;
  double SHx;
  double Dx;
  double Bx;
  double SVx;
  int loop_ub_tmp;
  int i;
  double tmp_data[1];
  double Ex;
  if (modes.useTurnSlip) {
    zeta1 = cos(atan(tirParams->PDXP1 * (tirParams->PDXP2 * incrVar.dfz + 1.0) *
                     cos(atan(tirParams->PDXP3 * postProInputs->kappa)) *
                     tirParams->UNLOADED_RADIUS * postProInputs->phi));
  } else {
    zeta1 = 1.0;
  }

  Cx = tirParams->PCX1 * tirParams->LCX;
  SHx = incrVar.dpi * incrVar.dpi;
  *Fx0 = (tirParams->PDX1 + tirParams->PDX2 * incrVar.dfz) * ((tirParams->PPX3 *
    incrVar.dpi + 1.0) + tirParams->PPX4 * SHx) * (1.0 - tirParams->PDX3 *
    (postProInputs->gamma * postProInputs->gamma)) * starVar.LMUX_star;
  if (postProInputs->Fz == 0.0) {
    *Fx0 = 0.0;
  }

  *mux = *Fx0;
  Dx = *Fx0 * postProInputs->Fz * zeta1;
  *Kxk = postProInputs->Fz * (tirParams->PKX1 + tirParams->PKX2 * incrVar.dfz) *
    exp(tirParams->PKX3 * incrVar.dfz) * ((tirParams->PPX1 * incrVar.dpi + 1.0)
    + tirParams->PPX2 * SHx) * tirParams->LKX;
  *Fx0 = Dx;
  if (Dx < 0.0) {
    *Fx0 = -1.0;
  } else if (Dx > 0.0) {
    *Fx0 = 1.0;
  } else {
    if (Dx == 0.0) {
      *Fx0 = 0.0;
    }
  }

  if (*Fx0 == 0.0) {
    *Fx0 = 1.0;
  }

  Bx = *Kxk / (Cx * Dx + 1.0E-6 * *Fx0);
  SHx = (tirParams->PHX1 + tirParams->PHX2 * incrVar.dfz) * tirParams->LHX;
  SVx = postProInputs->Fz * (tirParams->PVX1 + tirParams->PVX2 * incrVar.dfz) *
    tirParams->LVX * primeVar.LMUX_prime * zeta1;
  if (modes.isLowSpeed) {
    loop_ub_tmp = internalParams->reductionSmooth.size[0] *
      internalParams->reductionSmooth.size[1];
    for (i = 0; i < loop_ub_tmp; i++) {
      tmp_data[i] = SVx * internalParams->reductionSmooth.data[i];
    }

    *Fx0 = SVx;
    if (modes.isLowSpeed) {
      *Fx0 = tmp_data[0];
    }

    SVx = *Fx0;
    for (i = 0; i < loop_ub_tmp; i++) {
      tmp_data[i] = SHx * internalParams->reductionSmooth.data[i];
    }

    *Fx0 = SHx;
    if (modes.isLowSpeed) {
      *Fx0 = tmp_data[0];
    }

    SHx = *Fx0;
  }

  zeta1 = postProInputs->kappa + SHx;
  SHx = zeta1;
  if (zeta1 < 0.0) {
    SHx = -1.0;
  } else if (zeta1 > 0.0) {
    SHx = 1.0;
  } else {
    if (zeta1 == 0.0) {
      SHx = 0.0;
    }
  }

  Ex = ((tirParams->PEX1 + tirParams->PEX2 * incrVar.dfz) + tirParams->PEX3 *
        (incrVar.dfz * incrVar.dfz)) * (1.0 - tirParams->PEX4 * SHx) *
    tirParams->LEX;
  if (modes.useLimitsCheck && (Ex > 1.0)) {
    Ex = 1.0;
  }

  SHx = Bx * zeta1;
  *Fx0 = Dx * sin(Cx * atan(SHx - Ex * (SHx - atan(SHx)))) + SVx;
  loop_ub_tmp = 0;
  if (postProInputs->uVcx < 0.0) {
    loop_ub_tmp = 1;
  }

  if (0 <= loop_ub_tmp - 1) {
    tmp_data[0] = -*Fx0;
  }

  if (postProInputs->uVcx < 0.0) {
    *Fx0 = tmp_data[0];
  }
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                const i_struct_T *internalParams
 *                const b_struct_T modes
 *                const e_struct_T starVar
 *                const d_struct_T primeVar
 *                const c_struct_T incrVar
 *                double *Fy0
 *                double *muy
 *                double *Kya
 *                double *Kyg0
 *                double *SHy
 *                double *SVy
 *                double *By
 *                double *Cy
 * Return Type  : void
 */
static void Solver_calculateFy0(const struct0_T *tirParams, const struct_T
  *postProInputs, const i_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fy0, double *muy, double *Kya, double *Kyg0, double *SHy, double *SVy,
  double *By, double *Cy)
{
  double PHY3;
  double zeta2;
  double zeta3;
  double zeta4;
  double Dy;
  double SVyg;
  int zeta0;
  double signDy;
  int i;
  double tmp_data[1];
  double EHyp;
  PHY3 = 0.0;
  if (modes.useTurnSlip) {
    zeta3 = cos(atan(tirParams->PKYP1 * (tirParams->UNLOADED_RADIUS *
      tirParams->UNLOADED_RADIUS) * (postProInputs->phi * postProInputs->phi)));
    zeta4 = tirParams->UNLOADED_RADIUS * fabs(postProInputs->phi);
    zeta2 = cos(atan(tirParams->PDYP1 * (tirParams->PDYP2 * incrVar.dfz + 1.0) *
                     cos(atan(tirParams->PDYP3 * tan(postProInputs->alpha))) *
                     (zeta4 + tirParams->PDYP4 * sqrt(zeta4))));
  } else {
    zeta2 = 1.0;
    zeta3 = 1.0;
  }

  Dy = tirParams->PKY1 * primeVar.Fz0_prime;
  *Kya = Dy * (tirParams->PPY1 * incrVar.dpi + 1.0) * (1.0 - tirParams->PKY3 *
    0.0) * sin(tirParams->PKY4 * atan(postProInputs->Fz / primeVar.Fz0_prime /
    ((tirParams->PKY2 + tirParams->PKY5 * 0.0) * (tirParams->PPY2 * incrVar.dpi
    + 1.0)))) * zeta3 * tirParams->LKY;
  zeta4 = postProInputs->Fz * (tirParams->PVY3 + tirParams->PVY4 * incrVar.dfz);
  SVyg = zeta4 * 0.0 * tirParams->LKYC * primeVar.LMUY_prime * zeta2;
  if (tirParams->FITTYP == 6.0) {
    PHY3 = tirParams->PHY3;
    *Kyg0 = (tirParams->PHY3 * (Dy * sin(tirParams->PKY4 * atan
               (postProInputs->Fz / (tirParams->PKY2 * primeVar.Fz0_prime))) *
              tirParams->LKYC) + zeta4) * tirParams->LKYC;
  } else {
    *Kyg0 = postProInputs->Fz * (tirParams->PKY6 + tirParams->PKY7 * incrVar.dfz)
      * (tirParams->PPY5 * incrVar.dpi + 1.0) * tirParams->LKYC;
  }

  if (modes.useTurnSlip) {
    signDy = *Kya;
    if (*Kya < 0.0) {
      signDy = -1.0;
    } else if (*Kya > 0.0) {
      signDy = 1.0;
    } else {
      if (*Kya == 0.0) {
        signDy = 0.0;
      }
    }

    if (signDy == 0.0) {
      signDy = 1.0;
    }

    zeta3 = *Kya;
    if (*Kya < 0.0) {
      zeta3 = -1.0;
    } else if (*Kya > 0.0) {
      zeta3 = 1.0;
    } else {
      if (*Kya == 0.0) {
        zeta3 = 0.0;
      }
    }

    if (zeta3 == 0.0) {
      zeta3 = 1.0;
    }

    Dy = postProInputs->uVcx;
    if (postProInputs->uVcx < 0.0) {
      Dy = -1.0;
    } else if (postProInputs->uVcx > 0.0) {
      Dy = 1.0;
    } else {
      if (postProInputs->uVcx == 0.0) {
        Dy = 0.0;
      }
    }

    zeta4 = (tirParams->PHYP2 + tirParams->PHYP3 * incrVar.dfz) * Dy;
    EHyp = tirParams->PHYP4;
    if (modes.useLimitsCheck && (tirParams->PHYP4 > 1.0)) {
      EHyp = 1.0;
    }

    zeta0 = 0;
    zeta3 = *Kyg0 / (1.0 - internalParams->epsilong) / (tirParams->PHYP1 * zeta4
      * (*Kya + 1.0E-6 * zeta3)) * tirParams->UNLOADED_RADIUS *
      postProInputs->phi;
    zeta4 = (zeta4 * sin(tirParams->PHYP1 * atan(zeta3 - EHyp * (zeta3 - atan
                (zeta3)))) * Dy + 1.0) - SVyg / (*Kya + 1.0E-6 * signDy);
  } else {
    zeta0 = 1;
    zeta4 = 1.0;
  }

  signDy = *Kya;
  if (*Kya < 0.0) {
    signDy = -1.0;
  } else if (*Kya > 0.0) {
    signDy = 1.0;
  } else {
    if (*Kya == 0.0) {
      signDy = 0.0;
    }
  }

  if (signDy == 0.0) {
    signDy = 1.0;
  }

  if (tirParams->FITTYP == 6.0) {
    *SHy = (tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) * tirParams->LHY +
      PHY3 * 0.0 * tirParams->LKYC;
  } else {
    *SHy = (((tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) * tirParams->LHY
             + (*Kyg0 * 0.0 - SVyg) / (*Kya + 1.0E-6 * signDy) * (double)zeta0)
            + zeta4) - 1.0;
  }

  *SVy = postProInputs->Fz * (tirParams->PVY1 + tirParams->PVY2 * incrVar.dfz) *
    tirParams->LVY * primeVar.LMUY_prime * zeta2 + SVyg;
  if (modes.isLowSpeed) {
    zeta0 = internalParams->reductionSmooth.size[0] *
      internalParams->reductionSmooth.size[1];
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = *SVy * internalParams->reductionSmooth.data[i];
    }

    signDy = *SVy;
    if (modes.isLowSpeed) {
      signDy = tmp_data[0];
    }

    *SVy = signDy;
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = *SHy * internalParams->reductionSmooth.data[i];
    }

    signDy = *SHy;
    if (modes.isLowSpeed) {
      signDy = tmp_data[0];
    }

    *SHy = signDy;
  }

  zeta3 = starVar.alpha_star + *SHy;
  *Cy = tirParams->PCY1 * tirParams->LCY;
  *muy = (tirParams->PDY1 + tirParams->PDY2 * incrVar.dfz) * ((tirParams->PPY3 *
    incrVar.dpi + 1.0) + tirParams->PPY4 * (incrVar.dpi * incrVar.dpi)) * (1.0 -
    tirParams->PDY3 * 0.0) * starVar.LMUY_star;
  Dy = *muy * postProInputs->Fz * zeta2;
  signDy = zeta3;
  if (zeta3 < 0.0) {
    signDy = -1.0;
  } else if (zeta3 > 0.0) {
    signDy = 1.0;
  } else {
    if (zeta3 == 0.0) {
      signDy = 0.0;
    }
  }

  if (signDy == 0.0) {
    signDy = 1.0;
  }

  zeta4 = (tirParams->PEY1 + tirParams->PEY2 * incrVar.dfz) * ((tirParams->PEY5 *
    0.0 + 1.0) - (tirParams->PEY3 + tirParams->PEY4 * 0.0) * signDy) *
    tirParams->LEY;
  if (modes.useLimitsCheck && (zeta4 > 1.0)) {
    zeta4 = 1.0;
  }

  signDy = Dy;
  if (Dy < 0.0) {
    signDy = -1.0;
  } else if (Dy > 0.0) {
    signDy = 1.0;
  } else {
    if (Dy == 0.0) {
      signDy = 0.0;
    }
  }

  if (signDy == 0.0) {
    signDy = 1.0;
  }

  *By = *Kya / (*Cy * Dy + 1.0E-6 * signDy);
  zeta3 *= *By;
  *Fy0 = Dy * sin(*Cy * atan(zeta3 - zeta4 * (zeta3 - atan(zeta3)))) + *SVy;
  if (modes.useAlphaStar) {
    zeta0 = 0;
    if (postProInputs->uVcx < 0.0) {
      zeta0 = 1;
    }

    if (0 <= zeta0 - 1) {
      tmp_data[0] = -*Fy0;
    }

    signDy = *Fy0;
    if (postProInputs->uVcx < 0.0) {
      signDy = tmp_data[0];
    }

    *Fy0 = signDy;
  }

  signDy = *muy;
  if (postProInputs->Fz == 0.0) {
    signDy = 0.0;
  }

  *muy = signDy;
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                const h_struct_T *internalParams
 *                const b_struct_T modes
 *                const e_struct_T starVar
 *                const d_struct_T primeVar
 *                const c_struct_T incrVar
 *                double alphar
 *                double alphat
 *                double Kxk
 *                double Kya_prime
 *                double Fy
 *                double Fx
 *                double Dr
 *                double Cr
 *                double Br
 *                double Dt
 *                double Ct
 *                double Bt
 *                double Et
 *                double *Mz
 *                double *t
 *                double *Mzr
 * Return Type  : void
 */
static void Solver_calculateMz(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double alphar, double alphat, double Kxk, double Kya_prime, double Fy, double
  Fx, double Dr, double Cr, double Br, double Dt, double Ct, double Bt, double
  Et, double *Mz, double *t, double *Mzr)
{
  double a;
  double a_tmp;
  double SHyk;
  double b_a;
  double x;
  double s;
  double Mzr_tmp;
  int c_internalParams_sub0_reduction;
  int loop_ub_tmp;
  double d_internalParams_sub0_reduction[1];
  int e_internalParams_sub0_reduction;
  int loop_ub;
  double f_internalParams_sub0_reduction[1];
  int g_internalParams_sub0_reduction;
  double h_internalParams_sub0_reduction[1];
  int i_internalParams_sub0_reduction;
  double j_internalParams_sub0_reduction[1];
  struct_T postProInputs_sub0;
  e_struct_T starVar_sub0;
  i_struct_T expl_temp;
  double Fyo_sub0;
  double Byk;
  double Eyk;
  double unusedU2a;
  double unusedU2b;
  a = tan(alphar);
  a_tmp = Kxk / Kya_prime;
  SHyk = alphar;
  if (alphar < 0.0) {
    SHyk = -1.0;
  } else if (alphar > 0.0) {
    SHyk = 1.0;
  } else {
    if (alphar == 0.0) {
      SHyk = 0.0;
    }
  }

  b_a = tan(alphat);
  x = alphat;
  if (alphat < 0.0) {
    x = -1.0;
  } else if (alphat > 0.0) {
    x = 1.0;
  } else {
    if (alphat == 0.0) {
      x = 0.0;
    }
  }

  s = tirParams->UNLOADED_RADIUS * ((tirParams->SSZ1 + tirParams->SSZ2 * (Fy /
    tirParams->FNOMIN)) + (tirParams->SSZ3 + tirParams->SSZ4 * incrVar.dfz) *
    postProInputs->gamma) * tirParams->LS;
  Mzr_tmp = a_tmp * a_tmp * (postProInputs->kappa * postProInputs->kappa);
  *Mzr = Dr * cos(Cr * atan(Br * (atan(sqrt(a * a + Mzr_tmp)) * SHyk)));
  c_internalParams_sub0_reduction = internalParams->reductionSmooth.size[1];
  loop_ub_tmp = internalParams->reductionSmooth.size[0] *
    internalParams->reductionSmooth.size[1];
  if (0 <= loop_ub_tmp - 1) {
    memcpy(&d_internalParams_sub0_reduction[0],
           &internalParams->reductionSmooth.data[0], loop_ub_tmp * sizeof(double));
  }

  e_internalParams_sub0_reduction = internalParams->reductionSharp.size[1];
  loop_ub = internalParams->reductionSharp.size[0] *
    internalParams->reductionSharp.size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&f_internalParams_sub0_reduction[0],
           &internalParams->reductionSharp.data[0], loop_ub * sizeof(double));
  }

  g_internalParams_sub0_reduction = internalParams->reductionLinear.size[1];
  loop_ub = internalParams->reductionLinear.size[0] *
    internalParams->reductionLinear.size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&h_internalParams_sub0_reduction[0],
           &internalParams->reductionLinear.data[0], loop_ub * sizeof(double));
  }

  i_internalParams_sub0_reduction = internalParams->reductionLinear_alpha.size[1];
  loop_ub = internalParams->reductionLinear_alpha.size[0] *
    internalParams->reductionLinear_alpha.size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&j_internalParams_sub0_reduction[0],
           &internalParams->reductionLinear_alpha.data[0], loop_ub * sizeof
           (double));
  }

  postProInputs_sub0 = *postProInputs;
  postProInputs_sub0.gamma = 0.0;
  starVar_sub0 = starVar;
  starVar_sub0.gamma_star = 0.0;
  expl_temp.zeta8 = 1.0;
  expl_temp.zeta7 = 1.0;
  expl_temp.zeta6 = 1.0;
  expl_temp.zeta5 = 1.0;
  expl_temp.zeta4 = 1.0;
  expl_temp.zeta3 = 1.0;
  expl_temp.zeta1 = 1.0;
  expl_temp.zeta0 = 1.0;
  expl_temp.epsilong = internalParams->epsilong;
  expl_temp.zeta2 = 1.0;
  expl_temp.reductionLinear_alpha.size[0] = 1;
  expl_temp.reductionLinear_alpha.size[1] =
    internalParams->reductionLinear_alpha.size[1];
  if (0 <= i_internalParams_sub0_reduction - 1) {
    memcpy(&expl_temp.reductionLinear_alpha.data[0],
           &j_internalParams_sub0_reduction[0], i_internalParams_sub0_reduction *
           sizeof(double));
  }

  expl_temp.reductionLinear.size[0] = 1;
  expl_temp.reductionLinear.size[1] = internalParams->reductionLinear.size[1];
  if (0 <= g_internalParams_sub0_reduction - 1) {
    memcpy(&expl_temp.reductionLinear.data[0], &h_internalParams_sub0_reduction
           [0], g_internalParams_sub0_reduction * sizeof(double));
  }

  expl_temp.reductionSharp.size[0] = 1;
  expl_temp.reductionSharp.size[1] = internalParams->reductionSharp.size[1];
  if (0 <= e_internalParams_sub0_reduction - 1) {
    memcpy(&expl_temp.reductionSharp.data[0], &f_internalParams_sub0_reduction[0],
           e_internalParams_sub0_reduction * sizeof(double));
  }

  expl_temp.reductionSmooth.size[0] = 1;
  expl_temp.reductionSmooth.size[1] = internalParams->reductionSmooth.size[1];
  if (0 <= c_internalParams_sub0_reduction - 1) {
    memcpy(&expl_temp.reductionSmooth.data[0], &d_internalParams_sub0_reduction
           [0], c_internalParams_sub0_reduction * sizeof(double));
  }

  expl_temp.epsilonv = internalParams->epsilonv;
  expl_temp.epsilonr = internalParams->epsilonr;
  expl_temp.epsilony = internalParams->epsilony;
  expl_temp.epsilonk = internalParams->epsilonk;
  expl_temp.epsilonx = internalParams->epsilonx;
  Solver_calculateFy0(tirParams, &postProInputs_sub0, &expl_temp, modes,
                      starVar_sub0, primeVar, incrVar, &Fyo_sub0, &a_tmp, &a,
                      &SHyk, &Byk, &Eyk, &unusedU2a, &unusedU2b);
  a = a_tmp * postProInputs_sub0.Fz * ((tirParams->RVY1 + tirParams->RVY2 *
    incrVar.dfz) + tirParams->RVY3 * starVar.gamma_star) * cos(atan
    (tirParams->RVY4 * starVar.alpha_star)) * sin(tirParams->RVY5 * atan
    (tirParams->RVY6 * postProInputs_sub0.kappa)) * tirParams->LVYKA;
  SHyk = tirParams->RHY1 + tirParams->RHY2 * incrVar.dfz;
  Eyk = tirParams->REY1 + tirParams->REY2 * incrVar.dfz;
  if (modes.useLimitsCheck && (Eyk > 1.0)) {
    Eyk = 1.0;
  }

  Byk = (tirParams->RBY1 + tirParams->RBY4 * (starVar.gamma_star *
          starVar.gamma_star)) * cos(atan(tirParams->RBY2 * (starVar.alpha_star
    - tirParams->RBY3))) * tirParams->LYKA;
  if (modes.isLowSpeed) {
    for (c_internalParams_sub0_reduction = 0; c_internalParams_sub0_reduction <
         loop_ub_tmp; c_internalParams_sub0_reduction++) {
      d_internalParams_sub0_reduction[c_internalParams_sub0_reduction] = a *
        internalParams->reductionSmooth.data[c_internalParams_sub0_reduction];
    }

    a_tmp = a;
    if (modes.isLowSpeed) {
      a_tmp = d_internalParams_sub0_reduction[0];
    }

    a = a_tmp;
  }

  a_tmp = Bt * (atan(sqrt(b_a * b_a + Mzr_tmp)) * x);
  *t = Dt * cos(Ct * atan(a_tmp - Et * (a_tmp - atan(a_tmp)))) * cos
    (primeVar.alpha_prime) * tirParams->LFZO;
  if (modes.isLowSpeed) {
    for (c_internalParams_sub0_reduction = 0; c_internalParams_sub0_reduction <
         loop_ub_tmp; c_internalParams_sub0_reduction++) {
      d_internalParams_sub0_reduction[c_internalParams_sub0_reduction] = *t *
        internalParams->reductionSmooth.data[c_internalParams_sub0_reduction];
    }

    a_tmp = *t;
    if (modes.isLowSpeed) {
      a_tmp = d_internalParams_sub0_reduction[0];
    }

    *t = a_tmp;
    for (c_internalParams_sub0_reduction = 0; c_internalParams_sub0_reduction <
         loop_ub_tmp; c_internalParams_sub0_reduction++) {
      d_internalParams_sub0_reduction[c_internalParams_sub0_reduction] = *Mzr *
        internalParams->reductionSmooth.data[c_internalParams_sub0_reduction];
    }

    a_tmp = *Mzr;
    if (modes.isLowSpeed) {
      a_tmp = d_internalParams_sub0_reduction[0];
    }

    *Mzr = a_tmp;
  }

  if (tirParams->FITTYP == 6.0) {
    *Mz = (-*t * (Fy - a) + *Mzr) + s * Fx;
  } else {
    a = Byk * (postProInputs_sub0.kappa + SHyk);
    a_tmp = Byk * SHyk;
    *Mz = (-*t * (cos(tirParams->RCY1 * atan(a - Eyk * (a - atan(a)))) / cos
                  (tirParams->RCY1 * atan(a_tmp - Eyk * (a_tmp - atan(a_tmp)))) *
                  Fyo_sub0) + *Mzr) + s * Fx;
  }
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                const h_struct_T *internalParams
 *                const b_struct_T modes
 *                const e_struct_T starVar
 *                const d_struct_T primeVar
 *                const c_struct_T incrVar
 *                double Kya
 *                double SHy
 *                double SVy
 *                double By
 *                double Cy
 *                double *Mz0
 *                double *alphar
 *                double *alphat
 *                double *Dr
 *                double *Cr
 *                double *Br
 *                double *Dt
 *                double *Ct
 *                double *Bt
 *                double *Et
 *                double *Kya_prime
 * Return Type  : void
 */
static void Solver_calculateMz0(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  const e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double Kya, double SHy, double SVy, double By, double Cy, double *Mz0, double *
  alphar, double *alphat, double *Dr, double *Cr, double *Br, double *Dt, double
  *Ct, double *Bt, double *Et, double *Kya_prime)
{
  double x;
  double Fz;
  double Mzp_inf;
  double Bt_tmp;
  double alphay;
  double t0_tmp;
  double PHY3;
  double unusedU14;
  double SHyk;
  double SVyg;
  double Drp;
  int zeta0;
  int i;
  double tmp_data[1];
  double Fyo_sub0;
  double zeta2;
  double zeta6;
  h_struct_T unusedU20;
  double unusedU1f;
  x = postProInputs->uVcx;
  if (modes.useLimitsCheck) {
    Mzp_inf = postProInputs->Fz_lowLimit;
    if (postProInputs->Fz <= 0.0) {
      Mzp_inf = 0.0;
    }

    Fz = Mzp_inf;
  } else {
    Fz = postProInputs->Fz;
  }

  Mzp_inf = Kya;
  if (Kya < 0.0) {
    Mzp_inf = -1.0;
  } else if (Kya > 0.0) {
    Mzp_inf = 1.0;
  } else {
    if (Kya == 0.0) {
      Mzp_inf = 0.0;
    }
  }

  if (Mzp_inf == 0.0) {
    Mzp_inf = 1.0;
  }

  *Kya_prime = Kya + 1.0E-6 * Mzp_inf;
  *alphar = starVar.alpha_star + (SHy + SVy / *Kya_prime);
  *alphat = starVar.alpha_star + ((tirParams->QHZ1 + tirParams->QHZ2 *
    incrVar.dfz) + (tirParams->QHZ3 + tirParams->QHZ4 * incrVar.dfz) *
    starVar.gamma_star);
  if (modes.useTurnSlip) {
    Mzp_inf = cos(atan(tirParams->QDTP1 * tirParams->UNLOADED_RADIUS *
                       postProInputs->phi));
  } else {
    Mzp_inf = 1.0;
  }

  *Dt = (tirParams->QDZ1 + tirParams->QDZ2 * incrVar.dfz) * (1.0 -
    tirParams->PPZ1 * incrVar.dpi) * ((tirParams->QDZ3 * postProInputs->gamma +
    1.0) + tirParams->QDZ4 * (postProInputs->gamma * postProInputs->gamma)) * Fz
    * (tirParams->UNLOADED_RADIUS / primeVar.Fz0_prime) * tirParams->LTR *
    Mzp_inf;
  Bt_tmp = fabs(postProInputs->gamma);
  alphay = incrVar.dfz * incrVar.dfz;
  *Bt = ((tirParams->QBZ1 + tirParams->QBZ2 * incrVar.dfz) + tirParams->QBZ3 *
         alphay) * ((tirParams->QBZ4 * postProInputs->gamma + 1.0) +
                    tirParams->QBZ5 * Bt_tmp) * tirParams->LKY /
    starVar.LMUY_star;
  *Ct = tirParams->QCZ1;
  *Et = ((tirParams->QEZ1 + tirParams->QEZ2 * incrVar.dfz) + tirParams->QEZ3 *
         alphay) * ((tirParams->QEZ4 + tirParams->QEZ5 * starVar.gamma_star) *
                    0.63661977236758138 * atan(*Bt * tirParams->QCZ1 * *alphat)
                    + 1.0);
  if (modes.useLimitsCheck && (*Et > 1.0)) {
    *Et = 1.0;
  }

  t0_tmp = cos(primeVar.alpha_prime);
  PHY3 = 0.0;
  alphay = tirParams->PKY1 * primeVar.Fz0_prime;
  unusedU14 = alphay * (tirParams->PPY1 * incrVar.dpi + 1.0) * (1.0 -
    tirParams->PKY3 * 0.0) * sin(tirParams->PKY4 * atan(postProInputs->Fz /
    primeVar.Fz0_prime / ((tirParams->PKY2 + tirParams->PKY5 * 0.0) *
    (tirParams->PPY2 * incrVar.dpi + 1.0)))) * tirParams->LKY;
  SHyk = postProInputs->Fz * (tirParams->PVY3 + tirParams->PVY4 * incrVar.dfz);
  SVyg = SHyk * 0.0 * tirParams->LKYC * primeVar.LMUY_prime;
  if (tirParams->FITTYP == 6.0) {
    PHY3 = tirParams->PHY3;
    alphay = (tirParams->PHY3 * (alphay * sin(tirParams->PKY4 * atan
                (postProInputs->Fz / (tirParams->PKY2 * primeVar.Fz0_prime))) *
               tirParams->LKYC) + SHyk) * tirParams->LKYC;
  } else {
    alphay = postProInputs->Fz * (tirParams->PKY6 + tirParams->PKY7 *
      incrVar.dfz) * (tirParams->PPY5 * incrVar.dpi + 1.0) * tirParams->LKYC;
  }

  Mzp_inf = unusedU14;
  if (unusedU14 < 0.0) {
    Mzp_inf = -1.0;
  } else if (unusedU14 > 0.0) {
    Mzp_inf = 1.0;
  } else {
    if (unusedU14 == 0.0) {
      Mzp_inf = 0.0;
    }
  }

  if (Mzp_inf == 0.0) {
    Mzp_inf = 1.0;
  }

  if (tirParams->FITTYP == 6.0) {
    alphay = (tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) * tirParams->LHY
      + PHY3 * 0.0 * tirParams->LKYC;
  } else {
    alphay = (((tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) *
               tirParams->LHY + (alphay * 0.0 - SVyg) / (unusedU14 + 1.0E-6 *
                Mzp_inf)) + 1.0) - 1.0;
  }

  Drp = postProInputs->Fz * (tirParams->PVY1 + tirParams->PVY2 * incrVar.dfz) *
    tirParams->LVY * primeVar.LMUY_prime + SVyg;
  if (modes.isLowSpeed) {
    zeta0 = internalParams->reductionSmooth.size[0] *
      internalParams->reductionSmooth.size[1];
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = Drp * internalParams->reductionSmooth.data[i];
    }

    Mzp_inf = Drp;
    if (modes.isLowSpeed) {
      Mzp_inf = tmp_data[0];
    }

    Drp = Mzp_inf;
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = alphay * internalParams->reductionSmooth.data[i];
    }

    Mzp_inf = alphay;
    if (modes.isLowSpeed) {
      Mzp_inf = tmp_data[0];
    }

    alphay = Mzp_inf;
  }

  alphay += starVar.alpha_star;
  SHyk = tirParams->PCY1 * tirParams->LCY;
  PHY3 = (tirParams->PDY1 + tirParams->PDY2 * incrVar.dfz) * ((tirParams->PPY3 *
    incrVar.dpi + 1.0) + tirParams->PPY4 * (incrVar.dpi * incrVar.dpi)) * (1.0 -
    tirParams->PDY3 * 0.0) * starVar.LMUY_star * postProInputs->Fz;
  Mzp_inf = alphay;
  if (alphay < 0.0) {
    Mzp_inf = -1.0;
  } else if (alphay > 0.0) {
    Mzp_inf = 1.0;
  } else {
    if (alphay == 0.0) {
      Mzp_inf = 0.0;
    }
  }

  if (Mzp_inf == 0.0) {
    Mzp_inf = 1.0;
  }

  SVyg = (tirParams->PEY1 + tirParams->PEY2 * incrVar.dfz) * ((tirParams->PEY5 *
    0.0 + 1.0) - (tirParams->PEY3 + tirParams->PEY4 * 0.0) * Mzp_inf) *
    tirParams->LEY;
  if (modes.useLimitsCheck && (SVyg > 1.0)) {
    SVyg = 1.0;
  }

  Mzp_inf = PHY3;
  if (PHY3 < 0.0) {
    Mzp_inf = -1.0;
  } else if (PHY3 > 0.0) {
    Mzp_inf = 1.0;
  } else {
    if (PHY3 == 0.0) {
      Mzp_inf = 0.0;
    }
  }

  if (Mzp_inf == 0.0) {
    Mzp_inf = 1.0;
  }

  alphay *= unusedU14 / (SHyk * PHY3 + 1.0E-6 * Mzp_inf);
  Fyo_sub0 = PHY3 * sin(SHyk * atan(alphay - SVyg * (alphay - atan(alphay)))) +
    Drp;
  if (modes.useAlphaStar) {
    zeta0 = 0;
    if (postProInputs->uVcx < 0.0) {
      zeta0 = 1;
    }

    if (0 <= zeta0 - 1) {
      tmp_data[0] = -Fyo_sub0;
    }

    Mzp_inf = Fyo_sub0;
    if (postProInputs->uVcx < 0.0) {
      Mzp_inf = tmp_data[0];
    }

    Fyo_sub0 = Mzp_inf;
  }

  if (modes.useTurnSlip) {
    zeta0 = 0;
    zeta2 = internalParams->zeta2;
    zeta6 = cos(atan(tirParams->QBRP1 * tirParams->UNLOADED_RADIUS *
                     postProInputs->phi));
    unusedU20 = *internalParams;
    b_Solver_calculateFy0(tirParams, postProInputs, &unusedU20, modes, starVar,
                          primeVar, incrVar, &PHY3, &alphay, &SVyg, &SHyk, &Drp,
                          &unusedU14, &Mzp_inf, &unusedU1f);
    Mzp_inf = tirParams->QCRP1 * fabs(alphay) * tirParams->UNLOADED_RADIUS * Fz *
      sqrt(Fz / primeVar.Fz0_prime) * tirParams->LMP;
    if (Mzp_inf < 0.0) {
      Mzp_inf = 1.0E-6;
    }

    alphay = Mzp_inf / sin(1.5707963267948966 * tirParams->QDRP1);
    Drp = alphay * sin(tirParams->QDRP1 * atan(Fz * tirParams->UNLOADED_RADIUS *
      (tirParams->QDZ8 * tirParams->QDZ9 * incrVar.dfz + (tirParams->QDZ10 +
      tirParams->QDZ11 * incrVar.dfz * Bt_tmp)) * tirParams->LKZC /
      (tirParams->QDRP1 * alphay * (1.0 - internalParams->epsilong)) *
      tirParams->UNLOADED_RADIUS * postProInputs->phit));
    SHyk = tirParams->RHY1 + tirParams->RHY2 * incrVar.dfz;
    alphay = tirParams->REY1 + tirParams->REY2 * incrVar.dfz;
    if (modes.useLimitsCheck && (alphay > 1.0)) {
      alphay = 1.0;
    }

    PHY3 = (tirParams->RBY1 + tirParams->RBY4 * (starVar.gamma_star *
             starVar.gamma_star)) * cos(atan(tirParams->RBY2 *
      (starVar.alpha_star - tirParams->RBY3))) * tirParams->LYKA;
    SVyg = PHY3 * (postProInputs->kappa + SHyk);
    SHyk *= PHY3;
    *Cr = 0.63661977236758138 * acos(Mzp_inf * 0.63661977236758138 * atan
      (tirParams->QCRP2 * tirParams->UNLOADED_RADIUS * fabs(postProInputs->phit))
      * (cos(tirParams->RCY1 * atan(SVyg - alphay * (SVyg - atan(SVyg)))) / cos
         (tirParams->RCY1 * atan(SHyk - alphay * (SHyk - atan(SHyk))))) / fabs
      (Drp));
    alphay = Drp + 1.0;
  } else {
    zeta0 = 1;
    zeta2 = 1.0;
    zeta6 = 1.0;
    *Cr = 1.0;
    alphay = 1.0;
  }

  if (postProInputs->uVcx < 0.0) {
    x = -1.0;
  } else if (postProInputs->uVcx > 0.0) {
    x = 1.0;
  } else {
    if (postProInputs->uVcx == 0.0) {
      x = 0.0;
    }
  }

  *Dr = (Fz * tirParams->UNLOADED_RADIUS * ((tirParams->QDZ6 + tirParams->QDZ7 *
           incrVar.dfz) * tirParams->LRES * zeta2 + ((tirParams->QDZ8 +
            tirParams->QDZ9 * incrVar.dfz) * (tirParams->PPZ2 * incrVar.dpi +
            1.0) + (tirParams->QDZ10 + tirParams->QDZ11 * incrVar.dfz) * fabs
           (starVar.gamma_star)) * starVar.gamma_star * tirParams->LKZC *
          (double)zeta0) * starVar.LMUY_star * x * cos(starVar.alpha_star) +
         alphay) - 1.0;
  *Br = (tirParams->QBZ9 * tirParams->LKY / starVar.LMUY_star + tirParams->QBZ10
         * By * Cy) * zeta6;
  alphay = *Bt * *alphat;
  *Mz0 = -(*Dt * cos(tirParams->QCZ1 * atan(alphay - *Et * (alphay - atan(alphay))))
           * t0_tmp) * Fyo_sub0 + *Dr * cos(*Cr * atan(*Br * *alphar)) * t0_tmp;
}

/*
 * Arguments    : double tirParams_LONGVL
 *                double tirParams_UNLOADED_RADIUS
 *                double tirParams_FNOMIN
 *                double tirParams_VERTICAL_STIFFNESS
 *                double tirParams_BREFF
 *                double tirParams_DREFF
 *                double tirParams_FREFF
 *                double tirParams_Q_RE0
 *                double tirParams_Q_V1
 *                double tirParams_PFZ1
 *                double postProInputs_omega
 *                double postProInputs_uFz
 *                double postProInputs_ukappa
 *                double postProInputs_uVcx
 *                double dpi
 *                double *Re
 *                double *Romega
 *                double *omega
 * Return Type  : void
 */
static void Solver_calculateRe(double tirParams_LONGVL, double
  tirParams_UNLOADED_RADIUS, double tirParams_FNOMIN, double
  tirParams_VERTICAL_STIFFNESS, double tirParams_BREFF, double tirParams_DREFF,
  double tirParams_FREFF, double tirParams_Q_RE0, double tirParams_Q_V1, double
  tirParams_PFZ1, double postProInputs_omega, double postProInputs_uFz, double
  postProInputs_ukappa, double postProInputs_uVcx, double dpi, double *Re,
  double *Romega, double *omega)
{
  double a;
  double Cz;
  double Re_old;
  *omega = postProInputs_omega;
  a = postProInputs_omega * tirParams_UNLOADED_RADIUS / tirParams_LONGVL;
  *Romega = tirParams_UNLOADED_RADIUS * (tirParams_Q_RE0 + tirParams_Q_V1 * (a *
    a));
  *Re = tirParams_UNLOADED_RADIUS * 0.965;
  Cz = tirParams_VERTICAL_STIFFNESS * (tirParams_PFZ1 * dpi + 1.0);
  Re_old = tirParams_UNLOADED_RADIUS;
  while (fabs(Re_old - *Re) > 1.0E-9) {
    Re_old = *Re;
    *omega = (postProInputs_ukappa * postProInputs_uVcx + postProInputs_uVcx) / *
      Re;
    a = *omega * tirParams_UNLOADED_RADIUS / tirParams_LONGVL;
    *Romega = tirParams_UNLOADED_RADIUS * (tirParams_Q_RE0 + tirParams_Q_V1 * (a
      * a));
    a = postProInputs_uFz / tirParams_FNOMIN;
    *Re = *Romega - tirParams_FNOMIN / Cz * (tirParams_DREFF * atan
      (tirParams_BREFF * a) + tirParams_FREFF * a);
  }
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                h_struct_T *internalParams
 *                const b_struct_T modes
 *                const e_struct_T starVar
 *                const d_struct_T primeVar
 *                const c_struct_T incrVar
 *                double *Fy0
 *                double *muy
 *                double *Kya
 *                double *Kyg0
 *                double *SHy
 *                double *SVy
 *                double *By
 *                double *Cy
 * Return Type  : void
 */
static void b_Solver_calculateFy0(const struct0_T *tirParams, const struct_T
  *postProInputs, h_struct_T *internalParams, const b_struct_T modes, const
  e_struct_T starVar, const d_struct_T primeVar, const c_struct_T incrVar,
  double *Fy0, double *muy, double *Kya, double *Kyg0, double *SHy, double *SVy,
  double *By, double *Cy)
{
  double PHY3;
  double zeta2;
  double zeta3;
  double zeta4;
  double EHyp;
  double Dy;
  double DHyp;
  double Kya_tmp;
  double signKya0;
  double SVyg;
  int zeta0;
  int i;
  double tmp_data[1];
  PHY3 = 0.0;
  if (modes.useTurnSlip) {
    zeta3 = cos(atan(tirParams->PKYP1 * (tirParams->UNLOADED_RADIUS *
      tirParams->UNLOADED_RADIUS) * (postProInputs->phi * postProInputs->phi)));
    zeta4 = tirParams->UNLOADED_RADIUS * fabs(postProInputs->phi);
    zeta2 = cos(atan(tirParams->PDYP1 * (tirParams->PDYP2 * incrVar.dfz + 1.0) *
                     cos(atan(tirParams->PDYP3 * tan(postProInputs->alpha))) *
                     (zeta4 + tirParams->PDYP4 * sqrt(zeta4))));
  } else {
    zeta2 = 1.0;
    zeta3 = 1.0;
  }

  internalParams->zeta2 = zeta2;
  zeta4 = tirParams->PKY1 * primeVar.Fz0_prime;
  EHyp = zeta4 * (tirParams->PPY1 * incrVar.dpi + 1.0);
  Dy = tirParams->PPY2 * incrVar.dpi + 1.0;
  DHyp = postProInputs->Fz / primeVar.Fz0_prime;
  Kya_tmp = starVar.gamma_star * starVar.gamma_star;
  *Kya = EHyp * (1.0 - tirParams->PKY3 * fabs(starVar.gamma_star)) * sin
    (tirParams->PKY4 * atan(DHyp / ((tirParams->PKY2 + tirParams->PKY5 * Kya_tmp)
       * Dy))) * zeta3 * tirParams->LKY;
  signKya0 = postProInputs->Fz * (tirParams->PVY3 + tirParams->PVY4 *
    incrVar.dfz);
  SVyg = signKya0 * starVar.gamma_star * tirParams->LKYC * primeVar.LMUY_prime *
    zeta2;
  if (tirParams->FITTYP == 6.0) {
    PHY3 = tirParams->PHY3;
    *Kyg0 = (tirParams->PHY3 * (zeta4 * sin(tirParams->PKY4 * atan
               (postProInputs->Fz / (tirParams->PKY2 * primeVar.Fz0_prime))) *
              tirParams->LKYC) + signKya0) * tirParams->LKYC;
  } else {
    *Kyg0 = postProInputs->Fz * (tirParams->PKY6 + tirParams->PKY7 * incrVar.dfz)
      * (tirParams->PPY5 * incrVar.dpi + 1.0) * tirParams->LKYC;
  }

  if (modes.useTurnSlip) {
    zeta4 = EHyp * (1.0 - tirParams->PKY3 * 0.0) * sin(tirParams->PKY4 * atan
      (DHyp / ((tirParams->PKY2 + tirParams->PKY5 * 0.0) * Dy))) * zeta3 *
      tirParams->LKY;
    zeta3 = *Kya;
    if (*Kya < 0.0) {
      zeta3 = -1.0;
    } else if (*Kya > 0.0) {
      zeta3 = 1.0;
    } else {
      if (*Kya == 0.0) {
        zeta3 = 0.0;
      }
    }

    if (zeta3 == 0.0) {
      zeta3 = 1.0;
    }

    signKya0 = zeta4;
    if (zeta4 < 0.0) {
      signKya0 = -1.0;
    } else if (zeta4 > 0.0) {
      signKya0 = 1.0;
    } else {
      if (zeta4 == 0.0) {
        signKya0 = 0.0;
      }
    }

    if (signKya0 == 0.0) {
      signKya0 = 1.0;
    }

    Dy = postProInputs->uVcx;
    if (postProInputs->uVcx < 0.0) {
      Dy = -1.0;
    } else if (postProInputs->uVcx > 0.0) {
      Dy = 1.0;
    } else {
      if (postProInputs->uVcx == 0.0) {
        Dy = 0.0;
      }
    }

    DHyp = (tirParams->PHYP2 + tirParams->PHYP3 * incrVar.dfz) * Dy;
    EHyp = tirParams->PHYP4;
    if (modes.useLimitsCheck && (tirParams->PHYP4 > 1.0)) {
      EHyp = 1.0;
    }

    zeta0 = 0;
    zeta4 = *Kyg0 / (1.0 - internalParams->epsilong) / (tirParams->PHYP1 * DHyp *
      (zeta4 + 1.0E-6 * signKya0)) * tirParams->UNLOADED_RADIUS *
      postProInputs->phi;
    zeta4 = (DHyp * sin(tirParams->PHYP1 * atan(zeta4 - EHyp * (zeta4 - atan
                (zeta4)))) * Dy + 1.0) - SVyg / (*Kya + 1.0E-6 * zeta3);
  } else {
    zeta0 = 1;
    zeta4 = 1.0;
  }

  zeta3 = *Kya;
  if (*Kya < 0.0) {
    zeta3 = -1.0;
  } else if (*Kya > 0.0) {
    zeta3 = 1.0;
  } else {
    if (*Kya == 0.0) {
      zeta3 = 0.0;
    }
  }

  if (zeta3 == 0.0) {
    zeta3 = 1.0;
  }

  if (tirParams->FITTYP == 6.0) {
    *SHy = (tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) * tirParams->LHY +
      PHY3 * starVar.gamma_star * tirParams->LKYC;
  } else {
    *SHy = (((tirParams->PHY1 + tirParams->PHY2 * incrVar.dfz) * tirParams->LHY
             + (*Kyg0 * starVar.gamma_star - SVyg) / (*Kya + 1.0E-6 * zeta3) *
             (double)zeta0) + zeta4) - 1.0;
  }

  *SVy = postProInputs->Fz * (tirParams->PVY1 + tirParams->PVY2 * incrVar.dfz) *
    tirParams->LVY * primeVar.LMUY_prime * zeta2 + SVyg;
  if (modes.isLowSpeed) {
    zeta0 = internalParams->reductionSmooth.size[0] *
      internalParams->reductionSmooth.size[1];
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = *SVy * internalParams->reductionSmooth.data[i];
    }

    zeta3 = *SVy;
    if (modes.isLowSpeed) {
      zeta3 = tmp_data[0];
    }

    *SVy = zeta3;
    for (i = 0; i < zeta0; i++) {
      tmp_data[i] = *SHy * internalParams->reductionSmooth.data[i];
    }

    zeta3 = *SHy;
    if (modes.isLowSpeed) {
      zeta3 = tmp_data[0];
    }

    *SHy = zeta3;
  }

  zeta4 = starVar.alpha_star + *SHy;
  *Cy = tirParams->PCY1 * tirParams->LCY;
  *muy = (tirParams->PDY1 + tirParams->PDY2 * incrVar.dfz) * ((tirParams->PPY3 *
    incrVar.dpi + 1.0) + tirParams->PPY4 * (incrVar.dpi * incrVar.dpi)) * (1.0 -
    tirParams->PDY3 * Kya_tmp) * starVar.LMUY_star;
  Dy = *muy * postProInputs->Fz * zeta2;
  zeta3 = zeta4;
  if (zeta4 < 0.0) {
    zeta3 = -1.0;
  } else if (zeta4 > 0.0) {
    zeta3 = 1.0;
  } else {
    if (zeta4 == 0.0) {
      zeta3 = 0.0;
    }
  }

  if (zeta3 == 0.0) {
    zeta3 = 1.0;
  }

  EHyp = (tirParams->PEY1 + tirParams->PEY2 * incrVar.dfz) * ((tirParams->PEY5 *
    Kya_tmp + 1.0) - (tirParams->PEY3 + tirParams->PEY4 * starVar.gamma_star) *
    zeta3) * tirParams->LEY;
  if (modes.useLimitsCheck && (EHyp > 1.0)) {
    EHyp = 1.0;
  }

  zeta3 = Dy;
  if (Dy < 0.0) {
    zeta3 = -1.0;
  } else if (Dy > 0.0) {
    zeta3 = 1.0;
  } else {
    if (Dy == 0.0) {
      zeta3 = 0.0;
    }
  }

  if (zeta3 == 0.0) {
    zeta3 = 1.0;
  }

  *By = *Kya / (*Cy * Dy + 1.0E-6 * zeta3);
  zeta4 *= *By;
  *Fy0 = Dy * sin(*Cy * atan(zeta4 - EHyp * (zeta4 - atan(zeta4)))) + *SVy;
  if (modes.useAlphaStar) {
    zeta0 = 0;
    if (postProInputs->uVcx < 0.0) {
      zeta0 = 1;
    }

    if (0 <= zeta0 - 1) {
      tmp_data[0] = -*Fy0;
    }

    zeta3 = *Fy0;
    if (postProInputs->uVcx < 0.0) {
      zeta3 = tmp_data[0];
    }

    *Fy0 = zeta3;
  }

  zeta3 = *muy;
  if (postProInputs->Fz == 0.0) {
    zeta3 = 0.0;
  }

  *muy = zeta3;
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                const h_struct_T *internalParams
 *                const b_struct_T modes
 *                f_struct_T *forcesAndmoments
 *                g_struct_T *varinf
 * Return Type  : void
 */
void Solver_doForcesAndMoments(const struct0_T *tirParams, const struct_T
  *postProInputs, const h_struct_T *internalParams, const b_struct_T modes,
  f_struct_T *forcesAndmoments, g_struct_T *varinf)
{
  double Vsx_tmp_tmp;
  double Vsx;
  double Exa;
  double Vsy;
  double Vs_tmp;
  double Vc;
  double Fz0_prime;
  double incrVar_dfz;
  double incrVar_dpi;
  double starVar_alpha_star;
  double starVar_gamma_star;
  double signVc;
  double primeVar_alpha_prime;
  double LMUX_star;
  double LMUY_star;
  double primeVar_LMUX_prime;
  double primeVar_LMUY_prime;
  e_struct_T expl_temp;
  d_struct_T b_expl_temp;
  c_struct_T c_expl_temp;
  double Fx;
  double Kxk;
  h_struct_T b_internalParams;
  e_struct_T d_expl_temp;
  d_struct_T e_expl_temp;
  c_struct_T f_expl_temp;
  double Fy0;
  double muy;
  double Kya;
  e_struct_T g_expl_temp;
  d_struct_T h_expl_temp;
  c_struct_T i_expl_temp;
  double alphar;
  double alphat;
  double Dr;
  double Cr;
  double Br;
  double Dt;
  double Ct;
  double Bt;
  double Et;
  double Kya_prime;
  int loop_ub;
  int i;
  double x0_data[1];
  e_struct_T j_expl_temp;
  d_struct_T k_expl_temp;
  c_struct_T l_expl_temp;

  /* SOLVER Solver for Magic Formula 5.2, 6.1 and 6.2 Tyre Models */
  /*  All the functions available for the user are placed here. */
  /* FULLSTEADYSTATE Performs a full calculation of the MF model */
  /*  */
  /*  The calculation is split into 3 stages: */
  /*    1: The inputs are validated and analyzed by calling */
  /*    obj.parseInputs (public access) */
  /*  */
  /*    2: The forces and moments are calculated by calling */
  /*    obj.doForcesAndMoments (private access) */
  /*  */
  /*    3: Extra calculations like the lateral and longitudinal */
  /*    stiffness, contact patch length and width, etc by calling */
  /*    obj.doExtras (private access) */
  /*  */
  /*  Syntax: */
  /*  out = obj.fullSteadyState(tirParams, inputs, useMode) */
  /*  */
  /*  Returns a matrix with all the outputs for mfeval. */
  /*  For more information about the inputs, please refer to */
  /*  mfeval */
  /*  */
  /*  See also: mfeval */
  /*  Declare extrinsic functions for C code generation */
  /*  Solve in steady state mode */
  /*  Call parseInputs */
  /*  Call doForcesAndMoments */
  /*  Call doExtras */
  /*  Check the sign of the coefficient of friction */
  /*  The calculation of Fy is not affected by the sign of muy */
  /*  muy < 0 */
  /*  useLimitsCheck */
  /*  Preallocate out variable */
  /*  Pack all the outputs */
  /*  fullSteadyState */
  /* FULLCPI Performs the calculation of the MF model optimized */
  /* for the Simulink block with Contact Point Interface (CPI). */
  /* Less outputs are obtained compared to fullSteadyState */
  /*  */
  /*  The calculation is split into 3 sages: */
  /*    1: The inputs are validated and analyzed by calling */
  /*    obj.parseInputs (public access) */
  /*  */
  /*    2: The forces and moments are calculated by calling */
  /*    obj.doForcesAndMoments (private access) */
  /*  */
  /*    3: Extra calculations like the vertical deflection, loaded */
  /*    radius, etc are calculated by calling obj.doExtras (private */
  /*    access) */
  /*  */
  /*  Syntax: */
  /*  out = solver.obj(tirParams,inputs,userUseMode,userDynamics) */
  /*  */
  /*  Returns a matrix with the necessary outputs for the Simulink */
  /*  CPI block. */
  /*  */
  /*  Where: */
  /*    tirParams is a structure with the MF parameters */
  /*  */
  /*    inputs is an array of 9 elements: */
  /*    [Fz kappa alpha gamma phit Vx P omega Vsy] */
  /*  */
  /*    userUseMode is a number with values 1 or 2 to specify if */
  /*    tuenslip is applied (2) or not (1) */
  /*  */
  /*    userDynamics is a number between 1 to 3 to specify the */
  /*    dynamics [steady state (1), linear transients (2) or */
  /*    nonlinear (3)] */
  /*  */
  /*  See also: mfeval */
  /*  Get the right useMode number */
  /*  Use turn slip */
  /*  Do not use turn slip */
  /*  Parse Inputs */
  /*  Call doForcesAndMoments */
  /*  Call doExtras */
  /*  Preallocate OUTPUT variable */
  /*  Pack all the outputs */
  /*  fullCPI */
  /* PARSEINPUTS Applies certain limits on the inputs depending on */
  /* useMode and userDynamics to create necessary variables for */
  /* obj.fullSteadyState and obj.fullCPI */
  /*  */
  /*  Syntax: */
  /*  [postProInputs, internalParams, modes] = */
  /*  parseInputs(tirParams,inputs, useMode, userDynamics) */
  /*  */
  /*  Where: */
  /*  postProInputs is a structure with the original and */
  /*  post-processed inputs */
  /*  */
  /*  internalParams is a structure that contains parameters not */
  /*  defined in the TIR file with default values from the */
  /*  literature and internal variables for the low speed and */
  /*  turn slip modes. */
  /*  */
  /*  modes is a structure with logical variables that depend on */
  /*  the useMode number. */
  /*  */
  /*  See also: mfeval */
  /*  Declare extrinsic functions */
  /*  Check the number of digits of useMode */
  /*  Save the elements into separate variables */
  /*  Pre-allocate variables for C code generation */
  /*  Create the logical variables */
  /*  1: include limits checks */
  /*  2: revoke limits checks */
  /*  useLimitsCheck */
  /*  1: revoke alpha_star definition */
  /*  2: include alpha_star definition */
  /*  useAlphaStar */
  /*  1: combined force/moment calculation */
  /*  2: combined force/moment calculation + turn slip */
  /*  useAlphaStar */
  /*  Parameters not specified in the TIR file */
  /*  Used to avoid low speed singularity */
  /*  [Eqn (4.E6a) Page 178 - Book] */
  /*  Check the size of the inputs */
  /*  Flag an error if the number of inputs is incorrect */
  /*  if Dynamics is zero (MATLAB) */
  /*  Unpack the input matrix into separate input variables */
  /*  vertical load            (N) */
  /*  longitudinal slip        (-) (-1 = locked wheel) */
  /*  side slip angle          (radians) */
  /*  inclination angle        (radians) */
  /*  turn slip                (1/m) */
  /*  forward velocity         (m/s) */
  /*  IMPORTANT NOTE: Vx = Vcx [Eqn (7.4) Page 331 - Book] */
  /*  It is assumed that the difference between the wheel centre */
  /*  longitudinal velocity Vx and the longitudinal velocity Vcx of */
  /*  the contact centre is negligible */
  /*  If any Fz is negative set it to zero */
  /*  Create a copy of the variables (u stands for unlimited) */
  /*  If the pressure is specified by the user, grab it */
  /*  pressure (Pa) */
  /*  If is not, use the Inflation pressure of the TIR file */
  /*  pressure (Pa) */
  /*  if pressure */
  /*  Limits: */
  /*  This section applies the appropriate limits to the input */
  /*  values of the model based on the MF limits specified in */
  /*  the TIR File */
  /*  Fz_lowLimit is only used in some Moments equations */
  /*  Pre-declare the variable for C code generation */
  /*  Turn slip modifications */
  /*  Empirically discovered */
  /*  Minimum Speed */
  /*  Inflation Pressure Range */
  /*  Vertical Force Range */
  /*  Slip Angle Range */
  /*  Inclination Angle Range */
  /*  Long Slip Range */
  /*  Low Speed Model: */
  /*  Create a reduction factor for low speed and standstill */
  /*  Logic to flag if the speed it's below the limit */
  /*  Create a vector with numbers between 0 and 1 to apply a */
  /*  reduction factor with smooth transitions. */
  /*  Create a vector with numbers between 0 and 1 to apply a */
  /*  linear reduction toward zero */
  /*  Create a vector with numbers between 0 and 1 to apply a */
  /*  reduction factor using a sharp decrease toward zero but */
  /*  smooth transition toward VXLOW */
  /*  ukappaLow is equal to ukappa but with low speed */
  /*  correction. This is only used to export kappa */
  /*  If Vcx is close to zero, use linear reduction */
  /*  Calculate the lateral speed if userDynamics is zero */
  /*  (called from MATLAB) */
  /*  Dynamics is 1 or 3 (called from Simulink) */
  /*  Grab Vsy from inputs to avoid Vsy = 0 when Vx = 0 */
  /*  Vsy calculation */
  /*  If the speed is negative, the turn slip is also negative */
  /*  Sum the forward and lateral speeds */
  /*  The slip angle also suffers a reduction when the sum of */
  /*  Vx and Vy is less than VXLOW */
  /*  Create a vector with numbers between 0 and 1 to apply a */
  /*  linear reduction toward zero for alpha */
  /*  Solve only when is mfeval from MATLAB */
  /*  if userDynamics == 0 */
  /*  Solve for steady state (both MATLAB and Simulink) */
  /*  If Vcx is close to zero then SA is also 0, linear */
  /*  reduction */
  /*  if steady state */
  /*  Check Slip Angle limits */
  /*  Check camber limits */
  /*  Check Fz limits */
  /*  Create a copy of Fz and apply the low limit. */
  /*  This is only used in some Moments equations */
  /*  Check pressure limits */
  /*  Check slip ratio limits */
  /*  Flag if anything is out of range. */
  /*  Not using limits checks */
  /*  if useLimitsCheck */
  /*  Pack the outputs */
  /*  Set optional inputs to zero (needed for C code generation) */
  /*  Set optional inputs to zero (needed for C code generation) */
  /*  Zero speed (empirically discovered) */
  /*  [Eqn (4.E2b) Page 177 - Book] */
  /*  Check if omega is one input */
  /*  Grab omega if exists */
  /*  If omega is NOT one of the inputs: Estimate omega */
  /*  if omega */
  /*  Save omega */
  /*  Longitudinal slip velocity */
  /*  [Eqn (2.3) Page 64 - Book] */
  /*  [Eqn (4.E1) Page 177 - Book] */
  /*  [Eqn (4.E2a) Page 177 - Book] */
  /*  [Eqn (4.90) Page 186 - Book] Camber reduction factor */
  /*  Speed limits (avoid zero speed) */
  /* Vc_prime = Vcx; % From the book */
  /*  From the Equation Manual */
  /*  Singularity protected velocity, text on Page 184 */
  /*  Rearrange [Eqn (4.75) Page 183 - Book] */
  /*  [Eqn (4.76) Page 184 - Book] */
  /*  IMPORTANT NOTE: Eqn (4.76) has been modified */
  /*  In chapter 2.2 "Definition of tire input quantities" in the Pacejka */
  /*  book, it is assumed that the z-axis of the road coordinate system */
  /*  "points downward normal to the road plane" (p. 62). Due to this */
  /*  definition, Pacejka introduces the minus sign for the spin slip so */
  /*  that a positive spin leads to a positive torque Mz (p. 68).But in */
  /*  CarMaker (and other MBS software), the road coordinate system is */
  /*  orientated differently. The z-axis points upward to the */
  /*  road plane. Thus switching the signs is not needed here. */
  /*  if turnslip */
  /*  parseInputs */
  /*  Public methods */
  /*  This methods are public but they are hidden for the normal user */
  /*  Unpack Parameters */
  /* Nominal speed */
  /* Nominal tyre inflation pressure */
  /* Nominal wheel load */
  /*  Scale factor of nominal (rated) load */
  /*  Scale factor of Fx peak friction coefficient */
  /*  Scale factor of Fy peak friction coefficient */
  /*  New scaling factor in Pacejka 2012 with it's default value. */
  /*  This scaling factor is not present in the standard MF6.1 TIR */
  /*  files. */
  /*  Scale factor with slip speed Vs decaying friction */
  /*  Some basic calculations are done before calculating forces */
  /*  and moments */
  /*  Velocities in point S (slip point) */
  /*  [Eqn (4.E5) Page 181 - Book] */
  /*  [Eqn (2.12) Page 67 - Book] and [(4.E3) Page 177 - Book] */
  /*  Important Note: */
  /*  Due to the ISO sign convention, equation 2.12 does not need a */
  /*  negative sign. The Pacejka book is written in adapted SAE. */
  /*  [Eqn (3.39) Page 102 - Book] -> Slip velocity of the slip point S */
  /*  Velocities in point C (contact) */
  /*  Assumption from page 67 of the book, paragraph above Eqn (2.11) */
  /*  Velocity of the wheel contact centre C, Not described in the book but is the same as [Eqn (3.39) Page 102 - Book] */
  /*  Effect of having a tire with a different nominal load */
  /*  [Eqn (4.E1) Page 177 - Book] */
  /*  Normalized change in vertical load */
  /*  [Eqn (4.E2a) Page 177 - Book] */
  /*  Normalized change in inflation pressure */
  /*  [Eqn (4.E2b) Page 177 - Book] */
  /*  Use of star (*) definition. Only valid for the book */
  /*  implementation. TNO MF-Tyre does not use this. */
  /*  [Eqn (4.E3) Page 177 - Book] */
  /*  [Eqn (4.E4) Page 177 - Book] */
  /*  if useAlphaStar */
  /*  For the aligning torque at high slip angles */
  /*  [Eqn (4.E6a) Page 178 - Book] [sign(Vc) term explained on page 177] */
  /*  [Eqn (4.E6) Page 177 - Book] */
  /*  Slippery surface with friction decaying with increasing (slip) speed */
  /*  [Eqn (4.E7) Page 179 - Book] */
  /*  [Eqn (4.E7) Page 179 - Book] */
  /*  Digressive friction factor */
  /*  On Page 179 of the book is suggested Amu = 10, but after */
  /*  comparing the use of the scaling factors against TNO, Amu = 1 */
  /*  was giving perfect match */
  /*  [Eqn (4.E8) Page 179 - Book] */
  /*  [Eqn (4.E8) Page 179 - Book] */
  /*  Pack outputs */
  /*  calculateBasic */
  /*  Declare extrinsic functions */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of Fx shape factor */
  /*  Scale factor of Fx curvature factor */
  /*  Scale factor of Fx slip stiffness */
  /*  Scale factor of Fx horizontal shift */
  /*  Scale factor of Fx vertical shift */
  /* [LONGITUDINAL_COEFFICIENTS] */
  /* Shape factor Cfx for longitudinal force */
  /* Longitudinal friction Mux at Fznom */
  /* Variation of friction Mux with load */
  /* Variation of friction Mux with camber squared */
  /* Longitudinal curvature Efx at Fznom */
  /* Variation of curvature Efx with load */
  /* Variation of curvature Efx with load squared */
  /* Factor in curvature Efx while driving */
  /* Longitudinal slip stiffness Kfx./Fz at Fznom */
  /* Variation of slip stiffness Kfx./Fz with load */
  /* Exponent in slip stiffness Kfx./Fz with load */
  /* Horizontal shift Shx at Fznom */
  /* Variation of shift Shx with load */
  /* Vertical shift Svx./Fz at Fznom */
  /* Variation of shift Svx./Fz with load */
  /* linear influence of inflation pressure on longitudinal slip stiffness */
  /* quadratic influence of inflation pressure on longitudinal slip stiffness */
  /* linear influence of inflation pressure on peak longitudinal friction */
  /* quadratic influence of inflation pressure on peak longitudinal friction */
  /*  Unpack Parameters */
  /*  Turn slip */
  /*  [Eqn (4.106) Page 188 - Book] */
  /*  [Eqn (4.105) Page 188 - Book] */
  /*  if useTurnSlip */
  /*  (> 0) (4.E11) */
  /*  (4.E13) */
  /*  Zero Fz correction */
  /*  (> 0) (4.E12) */
  /*  (= BxCxDx = dFxo./dkx at kappax = 0) (= Cfk) (4.E15) */
  /*  If [Dx = 0] then [sign(0) = 0]. This is done to avoid [Kxk / 0 = NaN] in Eqn 4.E16 */
  /*  (4.E16) [sign(Dx) term explained on page 177] */
  /*  (4.E17) */
  /*  (4.E18) */
  /*  Low speed model */
  /*  if isLowSpeed */
  /*  (4.E10) */
  /*  Only in Linear Transients mode */
  /*  if linear transients */
  /*  (<=1) (4.E14) */
  /*  Limits check */
  /*  if Ex > 1 */
  /*  if useLimitsCheck */
  /*  Pure longitudinal force */
  /*  (4.E9) */
  /* if any(isLowSpeed) && (userDynamics==2) */
  /*     % Low speed model and Linear Transients */
  /*     Fx0(isLowSpeed) = Dx.*(1-internalParams.reductionSmooth)*sign(kappax)+ Fx0.*internalParams.reductionSmooth; */
  /* end % if linear transients and low speed */
  /*  Backward speed check */
  /*  if is not Nonlinear Transients */
  /*  calculateFx0 */
  /*  Declare extrinsic functions */
  /*  Unpack Parameters */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of Fy shape factor */
  /*  Scale factor of Fy curvature factor */
  /*  Scale factor of Fy cornering stiffness */
  /*  Scale factor of Fy horizontal shift */
  /*  Scale factor of Fy vertical shift */
  /*  Scale factor of camber force stiffness */
  /* [LATERAL_COEFFICIENTS] */
  /* Shape factor Cfy for lateral forces */
  /* Lateral friction Muy */
  /* Variation of friction Muy with load */
  /* Variation of friction Muy with squared camber */
  /* Lateral curvature Efy at Fznom */
  /* Variation of curvature Efy with load */
  /* Zero order camber dependency of curvature Efy */
  /* Variation of curvature Efy with camber */
  /* Variation of curvature Efy with camber squared */
  /* Maximum value of stiffness Kfy./Fznom */
  /* Load at which Kfy reaches maximum value */
  /* Variation of Kfy./Fznom with camber */
  /* Curvature of stiffness Kfy */
  /* Peak stiffness variation with camber squared */
  /* Fy camber stiffness factor */
  /* Vertical load dependency of camber stiffness */
  /* Horizontal shift Shy at Fznom */
  /* Variation of shift Shy with load */
  /* Vertical shift in Svy./Fz at Fznom */
  /* Variation of shift Svy./Fz with load */
  /* Variation of shift Svy./Fz with camber */
  /* Variation of shift Svy./Fz with camber and load */
  /* influence of inflation pressure on cornering stiffness */
  /* influence of inflation pressure on dependency of nominal tyre load on cornering stiffness */
  /* linear influence of inflation pressure on lateral peak friction */
  /* quadratic influence of inflation pressure on lateral peak friction */
  /* Influence of inflation pressure on camber stiffness */
  /*  Pre-allocate MF5.2 parameters for C code generation */
  /* Variation of shift Shy with camber */
  /*  Turn slip */
  /* [TURNSLIP_COEFFICIENTS] */
  /* Cornering stiffness reduction due to spin */
  /* Peak Fy reduction due to spin parameter */
  /* Peak Fy reduction due to spin with varying load parameter */
  /* Peak Fy reduction due to spin with alpha parameter */
  /* Peak Fy reduction due to square root of spin parameter */
  /* [DIMENSION] */
  /* Free tyre radius */
  /*  Inputs */
  /*  [Eqn (4.79) Page 185 - Book] */
  /*  [Eqn (4.78) Page 185 - Book] */
  /*  [Eqn (4.77) Page 184 - Book] */
  /*  No turn slip and small camber angles */
  /*  First paragraph on page 178 of the book */
  /*  if useTurnSlip calculate zeta2 and zeta3 */
  /*  Save zeta2 for Mz calculations */
  /*  (= ByCyDy = dFyo./dalphay at alphay = 0) (if gamma =0: =Kya0 = CFa) (PKY4=2)(4.E25) */
  /*  (4.E28) */
  /*  Check MF version */
  /*  Load MF5.2 parameters */
  /* Variation of shift Shy with camber */
  /*  MF5.2 equations from the MF-Tyre equation manual */
  /*  Simplified without pressure dependency */
  /*  MF6.1 and 6.2 equatons */
  /*  (=dFyo./dgamma at alpha = gamma = 0) (= CFgamma) (4.E30) */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of Fy cornering stiffness */
  /* [TURNSLIP_COEFFICIENTS] */
  /* Fy-alpha curve lateral shift limitation */
  /* Fy-alpha curve maximum lateral shift parameter */
  /* Fy-alpha curve maximum lateral shift varying with load parameter */
  /* Fy-alpha curve maximum lateral shift parameter */
  /* [LATERAL_COEFFICIENTS] */
  /* Maximum value of stiffness Kfy./Fznom */
  /* Load at which Kfy reaches maximum value */
  /* Variation of Kfy./Fznom with camber */
  /* Curvature of stiffness Kfy */
  /* Peak stiffness variation with camber squared */
  /* influence of inflation pressure on cornering stiffness */
  /* influence of inflation pressure on dependency of nominal tyre load on cornering stiffness */
  /* [DIMENSION] */
  /* Free tyre radius */
  /*  Inputs */
  /*  IMPORTANT NOTE: Explanation of the above equation, Kya0 */
  /*  Kya0 is the cornering stiffness when the camber angle is zero */
  /*  (gamma=0) which is again the product of the coefficients By, Cy and */
  /*  Dy at zero camber angle. Information from Kaustub Ragunathan, email: */
  /*  carmaker-service-uk@ipg-automotive.com */
  /*  If [Kya = 0] then [sign(0) = 0]. This is done to avoid [num / 0 = NaN] in Eqn 4.E39 */
  /*  If [Kya0 = 0] then [sign(0) = 0]. This is done to avoid [num / 0 = NaN] */
  /*  (4.E39) [sign(Kya) term explained on page 177] */
  /*  epsilonk is a small factor added to avoid the singularity condition during zero velocity (equation 308, CarMaker reference Manual). */
  /*  (>0) [Eqn (4.85) Page 186 - Book] */
  /*  [Eqn (4.86) Page 186 - Book] */
  /*  (<=1) [Eqn (4.87) Page 186 - Book] */
  /*  Limits check */
  /*  if EHyp > 1 */
  /*  if useLimitsCheck */
  /*  Eqn (4.89) */
  /* [Eqn (4.88) Page 186 - Book] */
  /*  [Eqn (4.80) Page 185 - Book] */
  /*  [Eqn (4.83) Page 186 - Book] */
  /*  [Eqn (4.84) Page 186 - Book] */
  /*  No turn slip and small camber angles */
  /*  First paragraph on page 178 of the book */
  /*  if useTurnSlip calculate zeta0 and zeta4 */
  /*  If [Kya = 0] then [sign(0) = 0]. This is done to avoid [num / 0 = NaN] in Eqn 4.E27 */
  /*  Check MF version */
  /*  MF5.2 equations */
  /*  From the MF-Tyre equation manual */
  /*  MF6.1 and 6.2 equatons */
  /*  (4.E27) [sign(Kya) term explained on page 177] */
  /*  (4.E29) */
  /*  Low speed model */
  /*  if isLowSpeed */
  /*  (4.E20) */
  /*  (> 0) (4.E21) */
  /*  (4.E23) */
  /*  (4.E22) */
  /*  (<=1)(4.E24) */
  /*  Limits check */
  /*  if Ey > 1 */
  /*  if useLimitsCheck */
  /*  If [Dy = 0] then [sign(0) = 0]. This is done to avoid [Kya / 0 = NaN] in Eqn 4.E26 */
  /*  (4.E26) [sign(Dy) term explained on page 177] */
  /*  (4.E19) */
  /*  Backward speed check for alpha_star */
  /*  if useAlphaStar */
  /* if any(isLowSpeed) && (modes.userDynamics==2) */
  /*     % Low speed model and Linear Transients */
  /*     Fy0(isLowSpeed) = Dy.*(1-internalParams.reductionSmooth).*sign(-alphay+eps)+ Fy0.*internalParams.reductionSmooth; */
  /* end % if linear transients and low speed */
  /*  Zero Fz correction */
  /*  calculateFy0 */
  /*  Declare extrinsic functions */
  /*  Unpack Parameters */
  /*  Set Fz to zero if the input is negative */
  /*  if useLimitsCheck */
  /* [DIMENSION] */
  /* Free tyre radius */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of Fy cornering stiffness */
  /*  Scale factor of peak of pneumatic trail */
  /*  Scale factor for offset of residual torque */
  /*  %Scale factor of camber torque stiffness */
  /* [ALIGNING_COEFFICIENTS] */
  /* Trail slope factor for trail Bpt at Fznom */
  /* Variation of slope Bpt with load */
  /* Variation of slope Bpt with load squared */
  /* Variation of slope Bpt with camber */
  /* Variation of slope Bpt with absolute camber */
  /* Factor for scaling factors of slope factor Br of Mzr */
  /* Factor for dimensionless cornering stiffness of Br of Mzr */
  /* Shape factor Cpt for pneumatic trail */
  /* Peak trail Dpt = Dpt.*(Fz./Fznom.*R0) */
  /* Variation of peak Dpt" with load */
  /* Variation of peak Dpt" with camber */
  /* Variation of peak Dpt" with camber squared */
  /* Peak residual torque Dmr" = Dmr./(Fz.*R0) */
  /* Variation of peak factor Dmr" with load */
  /* Variation of peak factor Dmr" with camber */
  /* Variation of peak factor Dmr" with camber and load */
  /* Variation of peak factor Dmr with camber squared */
  /* Variation of Dmr with camber squared and load */
  /* Trail curvature Ept at Fznom */
  /* Variation of curvature Ept with load */
  /* Variation of curvature Ept with load squared */
  /* Variation of curvature Ept with sign of Alpha-t */
  /* Variation of Ept with camber and sign Alpha-t */
  /* Trail horizontal shift Sht at Fznom */
  /* Variation of shift Sht with load */
  /* Variation of shift Sht with camber */
  /* Variation of shift Sht with camber and load */
  /* effect of inflation pressure on length of pneumatic trail */
  /* Influence of inflation pressure on residual aligning torque */
  /*  (4.E35) */
  /*  If [Kya = 0] then [sign(0) = 0]. This is done to avoid [num / 0 = NaN] in Eqn 4.E38 */
  /*  (4.E39) [sign(Kya) term explained on page 177] */
  /*  (4.E38) */
  /*  = alphaf (4.E37) */
  /*  (4.E34) */
  /* [Eqn (4.91) Page 186 - Book] */
  /*  if useTurnSlip */
  /*  Dt0 = Fz.*(R0./Fz0_prime).*(QDZ1 + QDZ2.*dfz).*(1 - PPZ1.*dpi).* LTR.*sign(Vcx); % (4.E42) */
  /*  Dt = Dt0.*(1 + QDZ3.*abs(gamma_star) + QDZ4.*gamma_star.^2).*zeta5; % (4.E43) */
  /*  */
  /*  IMPORTANT NOTE: The above original equation (4.E43) was not matching the */
  /*  TNO solver. The coefficient Dt affects the pneumatic trail (t) and the */
  /*  self aligning torque (Mz). */
  /*  It was observed that when negative inclination angles where used as */
  /*  inputs, there was a discrepancy between the TNO solver and mfeval. */
  /*  This difference comes from the term QDZ3, that in the original equation */
  /*  is multiplied by abs(gamma_star). But in the paper the equation is */
  /*  different and the abs() term is not written. Equation (A60) from the */
  /*  paper resulted into a perfect match with TNO. */
  /*  Keep in mind that the equations from the paper don't include turn slip */
  /*  effects. The term zeta5 has been added although it doesn't appear in the */
  /*  paper. */
  /*  Paper definition: */
  /*  (A60) */
  /*  Bt = (QBZ1 + QBZ2.*dfz + QBZ3.*dfz.^2).*(1 + QBZ5.*abs(gamma_star) + QBZ6.*gamma_star.^2).*LKY./LMUY_star; %(> 0)(4.E40) */
  /*  */
  /*  IMPORTANT NOTE: In the above original equation (4.E40) it is used the */
  /*  parameter QBZ6, which doesn't exist in the standard TIR files. Also note */
  /*  that on page 190 and 615 of the book a full set of parameters is given */
  /*  and QBZ6 doesn't appear. */
  /*  The equation has been replaced with equation (A58) from the paper. */
  /*  Paper definition: */
  /* (> 0) (A58) */
  /*  (> 0) (4.E41) */
  /*  (<=1) (4.E44) */
  /*  Limits check */
  /*  if Et > 1 */
  /*  if useLimitsCheck */
  /* t(aplhat)(4.E33) */
  /*  Evaluate Fy0 with gamma = 0 and phit = 0 */
  /*  gamma=phi=0 (4.E32) */
  /*  [Eqn (4.102) Page 188 - Book] */
  /*  [Eqn (4.95) Page 187 - Book] */
  /*  Mzp_inf should be always > 0 */
  /*  (>0) [Eqn (4.96) Page 187 - Book] */
  /*  [Eqn (4.94) Page 187 - Book] */
  /* [Eqn (4.99) Page 187 - Book] */
  /*  Eqn from the manual */
  /*  Eqn from the manual */
  /*  [Eqn (4.103) Page 188 - Book] */
  /*  Eqn from the manual */
  /*  if useTurnSlip */
  /*  (4.E47) */
  /*  preferred: qBz9 = 0 (4.E45) */
  /*  (4.E46) */
  /*  =Mzr(alphar)(4.E36) */
  /*  (4.E31) */
  /*  calculateMz0 */
  /*  Declare extrinsic functions */
  /*  Unpack Parameters */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of alpha influence on Fx */
  /* [LONGITUDINAL_COEFFICIENTS] */
  /* Slope factor for combined slip Fx reduction */
  /* Variation of slope Fx reduction with kappa */
  /* Influence of camber on stiffness for Fx combined */
  /* Shape factor for combined slip Fx reduction */
  /* Curvature factor of combined Fx */
  /* Curvature factor of combined Fx with load */
  /* Shift factor for combined slip Fx reduction */
  /*  (4.E55) */
  /*  (<= 1) (4.E56) */
  /*  Limits check */
  /*  if Exa > 1 */
  /*  if useLimitsCheck */
  /*  (4.E57) */
  /*  (> 0) (4.E54) */
  /*  (4.E53) */
  /*  (4.E52) */
  /*  (> 0)(4.E51 */
  /*  (4.E50) */
  /*  calculateFx */
  /*  Declare extrinsic functions */
  /*  if useTurnSlip */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of alpha influence on Fx */
  /*  Scale factor of kappa induced Fy */
  /* [LATERAL_COEFFICIENTS] */
  /* Slope factor for combined Fy reduction */
  /* Variation of slope Fy reduction with alpha */
  /* Shift term for alpha in slope Fy reduction */
  /* Influence of camber on stiffness of Fy combined */
  /* Shape factor for combined Fy reduction */
  /* Curvature factor of combined Fy */
  /* Curvature factor of combined Fy with load */
  /* Shift factor for combined Fy reduction */
  /* Shift factor for combined Fy reduction with load */
  /* Kappa induced side force Svyk./Muy.*Fz at Fznom */
  /* Variation of Svyk./Muy.*Fz with load */
  /* Variation of Svyk./Muy.*Fz with camber */
  /* Variation of Svyk./Muy.*Fz with alpha */
  /* Variation of Svyk./Muy.*Fz with kappa */
  /* Variation of Svyk./Muy.*Fz with atan(kappa) */
  /*  (4.E67) */
  /*  (4.E66) */
  /*  (4.E65) */
  /*  (<=1) (4.E64) */
  /*  Limits check */
  /*  if Eyk > 1 */
  /*  useLimitsCheck */
  /*  (4.E63) */
  /*  (> 0) (4.E62) */
  /*  (4.E61) */
  /*  (4.E60) */
  /*  (> 0)(4.E59) */
  /*  Low speed model */
  /*  Line for Simulink */
  /*  if isLowSpeed */
  /*  (4.E58) */
  /*  calculateFy */
  /*  Unpack Parameters */
  /*  Empirically discovered: */
  /*  If Fz is below FzMin a reduction factor is applied: */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* [DIMENSION] */
  /* Free tyre radius */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of Mx vertical shift */
  /*  Scale factor of overturning couple */
  /* [OVERTURNING_COEFFICIENTS] */
  /* Vertical shift of overturning moment */
  /* Camber induced overturning couple */
  /* Fy induced overturning couple */
  /* Mixed load lateral force and camber on Mx */
  /* Load effect on Mx with lateral force and camber */
  /* B-factor of load with Mx */
  /* Camber with load on Mx */
  /* Lateral force with load on Mx */
  /* B-factor of lateral force with load on Mx */
  /* Vertical force with camber on Mx */
  /* B-factor of vertical force with camber on Mx */
  /* Camber squared induced overturning moment */
  /* Lateral force induced overturning moment */
  /* Lateral force induced overturning moment with camber */
  /* Influence of inflation pressure on overturning moment */
  /*  Mx = R0.*Fz.*(QSX1.*LVMX - QSX2.*gamma.*(1 + PPMX1.*dpi) + QSX3.*((Fy)/Fz0)... */
  /*      + QSX4.*cos(QSX5.*atan((QSX6.*(Fz./Fz0)).^2)).*sin(QSX7.*gamma + QSX8.*atan... */
  /*      (QSX9.*((Fy)/Fz0))) + QSX10.*atan(QSX11.*(Fz./Fz0)).*gamma).*LMX; %(4.E69) */
  /*  */
  /*  IMPORTANT NOTE: The above original equation (4.E69) is not used */
  /*  because is not matching the results of the official TNO solver. Also, in */
  /*  the book definition and in the official paper, parameters QSX12 QSX13 and */
  /*  QSX14 are not used */
  /*  Instead of using equations (4.E69) from the book or (A47) from the */
  /*  official paper, it has been used the equation (49) described in the draft */
  /*  paper of Besselink (Not the official paper). This draft can be downloaded */
  /*  from: */
  /*  https://pure.tue.nl/ws/files/3139488/677330157969510.pdf */
  /*  purl.tue.nl/677330157969510.pdf */
  /*  Draft paper definition: */
  /*  (49) */
  /*  calculateMx */
  /*  Unpack Parameters */
  /*  Empirically discovered: */
  /*  If Fz is below FzMin a reduction factor is applied: */
  /* [OPERATING_CONDITIONS] */
  /* Nominal tyre inflation pressure */
  /* [MODEL] */
  /* Nominal speed */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* [DIMENSION] */
  /* Free tyre radius */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of rolling resistance torque */
  /* [ROLLING_COEFFICIENTS] */
  /* Rolling resistance torque coefficient */
  /* Rolling resistance torque depending on Fx */
  /* Rolling resistance torque depending on speed */
  /* Rolling resistance torque depending on speed .^4 */
  /* Rolling resistance torque depending on camber squared */
  /* Rolling resistance torque depending on load and camber squared */
  /* Rolling resistance torque coefficient load dependency */
  /* Rolling resistance torque coefficient pressure dependency */
  /*  My = Fz.R0*(QSY1 + QSY2.*(Fx./Fz0) + QSY3.*abs(Vcx./V0) + QSY4.*(Vcx./V0).^4 ... */
  /*      +(QSY5 + QSY6.*(Fz./Fz0)).*gamma.^2).*((Fz./Fz0).^QSY7.*(p./pi0).^QSY8).*LMY.; %(4.E70) */
  /*  */
  /*  IMPORTANT NOTE: The above equation from the book (4.E70) is not used */
  /*  because is not matching the results of the official TNO solver. */
  /*  This equation gives a positive output of rolling resistance, and in the */
  /*  ISO coordinate system, My should be negative. Furthermore, the equation */
  /*  from the book has an error, multiplying all the equation by Fz instead of */
  /*  Fz0 (first term). */
  /*  Because of the previous issues it has been used the equation (A48) from */
  /*  the paper. */
  /*  Check MF version */
  /*  MF5.2 equations */
  /*  From the MF-Tyre equation manual */
  /*  MF6.1 and MF6.2 equations */
  /*  Paper definition: */
  /* (A48) */
  /*  Backward speed check */
  /*  Low speed model (Empirically discovered) */
  /*  Points for the interpolation */
  /*  Call the interpolation function */
  /*  Reduce My values */
  /*  if idx */
  /*  Negative SR check */
  /*  Zero speed check */
  /*  calculateMy */
  /*  Unpack Parameters */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* [DIMENSION] */
  /* Free tyre radius */
  /* [SCALING_COEFFICIENTS] */
  /*  Scale factor of moment arm of Fx */
  /*  Scale factor of nominal (rated) load */
  /* [ALIGNING_COEFFICIENTS] */
  /* Nominal value of s./R0: effect of Fx on Mz */
  /* Variation of distance s./R0 with Fy./Fznom */
  /* Variation of distance s./R0 with camber */
  /* Variation of distance s./R0 with load and camber */
  /*  alphar_eq = sqrt(alphar.^2+(Kxk./Kya_prime).^2.*kappa.^2).*sign(alphar); % (4.E78) */
  /*  alphat_eq = sqrt(alphat.^2+(Kxk./Kya_prime).^2.*kappa.^2).*sign(alphat); % (4.E77) */
  /*  s = R0.*(SSZ1 + SSZ2.*(Fy./Fz0_prime) + (SSZ3 + SSZ4.*dfz).*gamma_star).*LS; % (4.E76) */
  /*  IMPORTANT NOTE: The equations 4.E78 and 4.E77 are not used due to small */
  /*  differences discovered at negative camber angles with the TNO solver. */
  /*  Instead equations A54 and A55 from the paper are used. */
  /*  */
  /*  IMPORTANT NOTE: The coefficient "s" (Equation 4.E76) determines the */
  /*  effect of Fx into Mz. The book uses "Fz0_prime" in the formulation, */
  /*  but the paper uses "Fz0". The equation (A56) from the paper has a better */
  /*  correlation with TNO. */
  /*  (A54) */
  /*  (A55) */
  /*  (A56) */
  /*  (4.E75) */
  /*  Evaluate Fy and Fy0 with gamma = 0 and phit = 0 */
  /*  Evaluate Fy0 with gamma = 0 and phit  = 0 */
  /*  Evaluate Gyk with gamma = 0 and phit  = 0 */
  /*  Note: in the above equation starVar is used instead of */
  /*  starVar_sub0 beacuse it was found a better match with TNO */
  /*  (4.E74) */
  /*  (4.E73) */
  /*  IMPORTANT NOTE: the above equation is not written in any source, but "t" */
  /*  is multiplied by LFZO in the TNO dteval function. This has been empirically */
  /*  discovered. */
  /*  Low speed model */
  /*  Line for Simulink */
  /*  if isLowSpeed */
  /*  (4.E72) */
  /*  Check MF version */
  /*  MF5.2 equations */
  /*  % From the MF-Tyre equation manual */
  /*  MF6.1 and 6.2 equatons */
  /*  (4.E71) */
  /*  calculateMz */
  /*  Unpack Parameters */
  /*  Rename the TIR file variables in the Pacejka style */
  /*  Nominal (rated) wheel load */
  /*  Free tyre radius */
  /*  Nominal speed (LONGVL) */
  /*  Vertical stiffness */
  /* Pressure effect on vertical stiffness */
  /* Low load stiffness effective rolling radius */
  /* Peak value of effective rolling radius */
  /* High load stiffness effective rolling radius */
  /* Ratio of free tyre radius with nominal tyre radius */
  /* Tyre radius increase with speed */
  /*  C code declaration */
  /*  rotational speed (rad/s) */
  /*  [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius */
  /*  Excerpt from OpenTIRE MF6.1 implementation */
  /*  Date: 2016-12-01 */
  /*  Prepared for Marco Furlan/JLR */
  /*  Questions: henning.olsson@calspan.com */
  /*  Nominal stiffness (pressure corrected) */
  /*  [Eqn (5) Page 2 - Paper] - Vertical stiffness adapted for tyre inflation pressure */
  /*  Check if omega is one of the inputs */
  /*  If it is, use it to calculate Re, otherwise it can be approximated with a */
  /*  short iteration */
  /*  Omega is one of the inputs */
  /*  rotational speed (rad/s) [eps is added to avoid Romega = 0] */
  /*  [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius */
  /*  Eff. Roll. Radius */
  /*  [Eqn (7) Page 2 - Paper] */
  /*  Omega is not specified and is going to be approximated */
  /*  Initial guess of Re based on something slightly less than R0 */
  /*  Use the most up to date Re to calculate an omega */
  /*  omega = Vcx ./ Re; % Old definition of Henning without kappa, not valid for brake and drive */
  /*  [Eqn (2.5) Page 65 - Book] */
  /*  Then we calculate free-spinning radius */
  /*  [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius */
  /*  Effective Rolling Radius */
  /*  [Eqn (7) Page 2 - Paper] */
  /*  while Re has not converged */
  /*  if omega is an input */
  /*  calculateRe */
  /*  Unpack Parameters */
  /*  Declare extrinsic functions */
  /* [MODEL] */
  /* Nominal speed */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* Tyre vertical stiffness */
  /* Vertical stiffness increase with speed */
  /* Quadratic term in load vs. deflection */
  /* Longitudinal force influence on vertical stiffness */
  /* Lateral force influence on vertical stiffness */
  /* Pressure effect on vertical stiffness */
  /* [DIMENSION] */
  /* Free tyre radius */
  /*  Model parameters as QFZ1 that normally aren't present in the TIR files */
  /*  Rearranging [Eqn (4) Page 2 - Paper] */
  /*  Split Eqn (A3.3) Page 619 of the Book into different bits: */
  /*  Joining all the effects except tyre deflection terms: */
  /*  Equation (A3.3) can be written as: */
  /*  Fz = (Q_FZ2*(rho/R0)^2 + Q_FZ1*(rho/R0)) * external_effects */
  /*  Rearranging all the terms we end up with a quadratic equation as: */
  /*  ax^2 + bx + c = 0 */
  /*  Q_FZ2*(rho/R0)^2 + Q_FZ1*(rho/R0) -(Fz/(external_effects)) = 0 */
  /*  Note: use of capital letters to avoid confusion with contact patch */
  /*  lengths "a" and "b" */
  /*  if there is a real solution */
  /*  tyre deflection for a free rolling tyre */
  /*  The loaded radius is the free-spinning radius minus the deflection */
  /*  Eqn A3.2 Page 619 - Book assuming positive rho at all the time */
  /*  Avoid division between zero */
  /*  Vertical stiffness (Spring) */
  /*  calculateRhoRl61 */
  /*  Pre-allocate the variable */
  /*  The loaded radius (Rl) cannot be calculated straight forward */
  /*  when the vertical load (Fz) is an input. */
  /*  */
  /*  Here I have written the equation for the vertical load (Fz) */
  /*  of the MF-Tyre 6.2 model that has the loaded radius (Rl) as */
  /*  input. (see calculateFz62 function) */
  /*  */
  /*  The Rl is calculated with a derivative-free method by */
  /*  minimizing the error between the target Fz (input to mfeval) */
  /*  and the Fz from  the "calculateFz62" function */
  /*  Declare the anonymous function (Cost function) for the fitting. */
  /*  The cost function will measure the error between your model and the data */
  /*  The @ operator creates the handle, and the parentheses () immediately */
  /*  after the @ operator include the function input arguments */
  /*  Set options for Newton-Raphson and Solve */
  /*  Call calculateFz62 with the calculated Rl to get rho_z */
  /*  Avoid division between zero */
  /*  Vertical stiffness (Spring) */
  /*  calculateRhoRl62 */
  /*  Unpack Parameters */
  /*  Rename the TIR file variables in the Pacejka style */
  /*  Free tyre radius */
  /*  Nominal width of the tyre */
  /*  Vertical stiffness */
  /* Pressure effect on vertical stiffness */
  /*  Nominal stiffness (pressure corrected) */
  /*  [Eqn (5) Page 2 - Paper] - Vertical stiffness adapted for tyre inflation pressure */
  /* [CONTACT_PATCH] */
  /* Square root term in contact length equation */
  /* Linear term in contact length equation */
  /* Root term in contact width equation */
  /* Linear term in contact width equation */
  /*  Pre-allocate MF5.2 parameters for C code generation */
  /* [DIMENSION] */
  /* Nominal rim radius */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* Distance to rim when bottoming starts to occur */
  /*  Approximated loaded Radius */
  /*  Bottoming model (Empirically discovered): */
  /*  Check if bottoming has happened */
  /*  Calculate the max Fz if bottoming happens to calculate the */
  /*  contact patch */
  /*  Substitute max Fz for the calculations */
  /*  Check MF version */
  /*  Load MF5.2 parameters */
  /*  MF5.2 Square root load term in contact length */
  /*  MF5.2 Linear load term in contact length */
  /*  Set default vaules (Empirically discovered) */
  /*  From the MF-Tyre equation manual */
  /*  From the MF-Tyre equation manual */
  /*  MF6.1 and 6.2 equatons */
  /* [Eqn (9) Page 3 - Paper] */
  /* [Eqn (10) Page 3 - Paper] */
  /*  calculateContactPatch */
  /*  Unpack Parameters */
  /* Tyre overall longitudinal stiffness vertical deflection dependency linear term */
  /* Tyre overall longitudinal stiffness vertical deflection dependency quadratic term */
  /* Tyre overall longitudinal stiffness pressure dependency */
  /* Tyre overall lateral stiffness vertical deflection dependency linear term */
  /* Tyre overall lateral stiffness vertical deflection dependency quadratic term */
  /* Tyre overall lateral stiffness pressure dependency */
  /*  Scale factor of nominal (rated) load */
  /*  Nominal (rated) wheel load */
  /*  Free tyre radius */
  /*  Reference pressure */
  /* Tyre overall longitudinal stiffness */
  /* Tyre overall lateral stiffness */
  /*  Pre-allocate MF5.2 parameters for C code generation */
  /*  Basic calculations */
  /*  [Eqn (4.E1) Page 177 - Book] */
  /*  [Eqn (4.E2a) Page 177 - Book] */
  /*  [Eqn (4.E2b) Page 177 - Book] */
  /*  Overall longitudinal Cx and lateral stiffness Cy */
  /*  (Eqn 17 - Paper) */
  /*  (Eqn 18 - Paper) */
  /*  Check MF version */
  /*  Load MF5.2 parameters */
  /*  Relaxation length SigKap0/Fz at Fznom */
  /*  Variation of SigKap0/Fz with load */
  /*  Variation of SigKap0/Fz with exponent of load */
  /*  Peak value of relaxation length SigAlp0/R0 */
  /*  Value of Fz/Fznom where SigAlp0 is extreme */
  /*  Scale factor of Relaxation length of Fx */
  /*  Scale factor of Relaxation length of Fy */
  /*  Variation of Kfy/Fznom with camber */
  /*  MF 5.2 equation for the longitudinal relaxation length */
  /*  From the MF-Tyre equation manual */
  /*  MF 5.2 equations for the lateral relaxation length */
  /*  From the MF-Tyre equation manual */
  /*  MF 6.1 and 6.2 equations for the relaxation lengths */
  /* (Eqn 19 - Paper) */
  /* (Eqn 19 - Paper) */
  /*  calculateRelax */
  /*  Unpack Parameters */
  /*  Derivative and append to get same number of elements */
  /*  calculateInstantaneousKya */
  /*  Calculate the vertical Force using the equations described in */
  /*  the MF-Tyre/MF-Swift 6.2 equation manual Document revision: */
  /*  20130706 */
  /* [MODEL] */
  /*  Nominal speed (LONGVL) */
  /* [DIMENSION] */
  /* Free tyre radius */
  /* [VERTICAL] */
  /* Nominal wheel load */
  /* Tyre vertical stiffness */
  /* Vertical stiffness increase with speed */
  /* Quadratic term in load vs. deflection */
  /* Longitudinal force influence on vertical stiffness */
  /* Lateral force influence on vertical stiffness */
  /* Pressure effect on vertical stiffness */
  /* Explicit load dependency for including the lateral force influence on vertical stiffness */
  /* Linear load dependent camber angle influence on vertical stiffness */
  /* Quadratic load dependent camber angle influence on vertical stiffness */
  /* Linear load and camber angle dependent reduction on vertical stiffness */
  /* Combined camber angle and side slip angle effect on vertical stiffness (constant) */
  /* Combined camber angle and side slip angle linear effect on vertical stiffness */
  /* Combined camber angle and side slip angle quadratic effect on vertical stiffness */
  /*  Model parameters as QFZ1 that normally aren't present in the TIR files */
  /*  Rearranging [Eqn (4) Page 2 - Paper] */
  /*  Asymmetric effect for combinations of camber and lateral force */
  /*  Tyre deflection for a free rolling tyre */
  /*  Reference tread width */
  /*  Deflection caused by camber */
  /*  Change NaN to Zero */
  /*  Vertical deflection */
  /*  Correction term */
  /*  Vertical force */
  /*  calculateFz62 */
  /*  Hidden methods */
  /*  This methods are private and only accessible inside the class */
  Vsx_tmp_tmp = fabs(postProInputs->uVcx);
  Vsx = -postProInputs->kappa * Vsx_tmp_tmp;
  Exa = tan(postProInputs->alpha);
  Vsy = Exa * Vsx_tmp_tmp;
  Vs_tmp = Vsy * Vsy;
  Vc = sqrt(postProInputs->uVcx * postProInputs->uVcx + Vs_tmp);
  Fz0_prime = tirParams->LFZO * tirParams->FNOMIN;
  incrVar_dfz = (postProInputs->Fz - Fz0_prime) / Fz0_prime;
  incrVar_dpi = (postProInputs->p - tirParams->NOMPRES) / tirParams->NOMPRES;
  if (modes.useAlphaStar) {
    Vsy = postProInputs->uVcx;
    if (postProInputs->uVcx < 0.0) {
      Vsy = -1.0;
    } else if (postProInputs->uVcx > 0.0) {
      Vsy = 1.0;
    } else {
      if (postProInputs->uVcx == 0.0) {
        Vsy = 0.0;
      }
    }

    starVar_alpha_star = Exa * Vsy;
    starVar_gamma_star = sin(postProInputs->gamma);
  } else {
    starVar_alpha_star = postProInputs->alpha;
    starVar_gamma_star = postProInputs->gamma;
  }

  signVc = Vc;
  if (Vc < 0.0) {
    signVc = -1.0;
  } else if (Vc > 0.0) {
    signVc = 1.0;
  } else {
    if (Vc == 0.0) {
      signVc = 0.0;
    }
  }

  if (signVc == 0.0) {
    signVc = 1.0;
  }

  primeVar_alpha_prime = acos(postProInputs->uVcx / (Vc + 1.0E-6 * signVc));
  Vsy = 0.0 * sqrt(Vsx * Vsx + Vs_tmp) / tirParams->LONGVL + 1.0;
  LMUX_star = tirParams->LMUX / Vsy;
  LMUY_star = tirParams->LMUY / Vsy;
  primeVar_LMUX_prime = LMUX_star / (0.0 * LMUX_star + 1.0);
  primeVar_LMUY_prime = LMUY_star / (0.0 * LMUY_star + 1.0);
  expl_temp.LMUY_star = LMUY_star;
  expl_temp.LMUX_star = LMUX_star;
  expl_temp.gamma_star = starVar_gamma_star;
  expl_temp.alpha_star = starVar_alpha_star;
  b_expl_temp.LMUY_prime = primeVar_LMUY_prime;
  b_expl_temp.LMUX_prime = primeVar_LMUX_prime;
  b_expl_temp.alpha_prime = primeVar_alpha_prime;
  b_expl_temp.Fz0_prime = Fz0_prime;
  c_expl_temp.dpi = incrVar_dpi;
  c_expl_temp.dfz = incrVar_dfz;
  Solver_calculateFx0(tirParams, postProInputs, internalParams, modes, expl_temp,
                      b_expl_temp, c_expl_temp, &Fx, &varinf->mux, &Kxk);
  b_internalParams = *internalParams;
  d_expl_temp.LMUY_star = LMUY_star;
  d_expl_temp.LMUX_star = LMUX_star;
  d_expl_temp.gamma_star = starVar_gamma_star;
  d_expl_temp.alpha_star = starVar_alpha_star;
  e_expl_temp.LMUY_prime = primeVar_LMUY_prime;
  e_expl_temp.LMUX_prime = primeVar_LMUX_prime;
  e_expl_temp.alpha_prime = primeVar_alpha_prime;
  e_expl_temp.Fz0_prime = Fz0_prime;
  f_expl_temp.dpi = incrVar_dpi;
  f_expl_temp.dfz = incrVar_dfz;
  b_Solver_calculateFy0(tirParams, postProInputs, &b_internalParams, modes,
                        d_expl_temp, e_expl_temp, f_expl_temp, &Fy0, &muy, &Kya,
                        &Vsy, &Vc, &Vs_tmp, &Vsx, &Exa);
  g_expl_temp.LMUY_star = LMUY_star;
  g_expl_temp.LMUX_star = LMUX_star;
  g_expl_temp.gamma_star = starVar_gamma_star;
  g_expl_temp.alpha_star = starVar_alpha_star;
  h_expl_temp.LMUY_prime = primeVar_LMUY_prime;
  h_expl_temp.LMUX_prime = primeVar_LMUX_prime;
  h_expl_temp.alpha_prime = primeVar_alpha_prime;
  h_expl_temp.Fz0_prime = Fz0_prime;
  i_expl_temp.dpi = incrVar_dpi;
  i_expl_temp.dfz = incrVar_dfz;
  Solver_calculateMz0(tirParams, postProInputs, &b_internalParams, modes,
                      g_expl_temp, h_expl_temp, i_expl_temp, Kya, Vc, Vs_tmp,
                      Vsx, Exa, &Vsy, &alphar, &alphat, &Dr, &Cr, &Br, &Dt, &Ct,
                      &Bt, &Et, &Kya_prime);
  Exa = tirParams->REX1 + tirParams->REX2 * incrVar_dfz;
  if (modes.useLimitsCheck && (Exa > 1.0)) {
    Exa = 1.0;
  }

  Vs_tmp = starVar_gamma_star * starVar_gamma_star;
  Vsy = (tirParams->RBX1 + tirParams->RBX3 * Vs_tmp) * cos(atan(tirParams->RBX2 *
    postProInputs->kappa)) * tirParams->LXAL;
  Vc = Vsy * (starVar_alpha_star + tirParams->RHX1);
  Vsy *= tirParams->RHX1;
  Fx *= cos(tirParams->RCX1 * atan(Vc - Exa * (Vc - atan(Vc)))) / cos
    (tirParams->RCX1 * atan(Vsy - Exa * (Vsy - atan(Vsy))));
  if (modes.useTurnSlip) {
    Vsx = b_internalParams.zeta2;
  } else {
    Vsx = 1.0;
  }

  signVc = muy * postProInputs->Fz * ((tirParams->RVY1 + tirParams->RVY2 *
    incrVar_dfz) + tirParams->RVY3 * starVar_gamma_star) * cos(atan
    (tirParams->RVY4 * starVar_alpha_star)) * Vsx * sin(tirParams->RVY5 * atan
    (tirParams->RVY6 * postProInputs->kappa)) * tirParams->LVYKA;
  Vc = tirParams->RHY1 + tirParams->RHY2 * incrVar_dfz;
  Vsx = tirParams->REY1 + tirParams->REY2 * incrVar_dfz;
  if (modes.useLimitsCheck && (Vsx > 1.0)) {
    Vsx = 1.0;
  }

  Vsy = (tirParams->RBY1 + tirParams->RBY4 * Vs_tmp) * cos(atan(tirParams->RBY2 *
    (starVar_alpha_star - tirParams->RBY3))) * tirParams->LYKA;
  if (modes.isLowSpeed) {
    loop_ub = b_internalParams.reductionSmooth.size[0] *
      b_internalParams.reductionSmooth.size[1];
    for (i = 0; i < loop_ub; i++) {
      x0_data[i] = signVc * b_internalParams.reductionSmooth.data[i];
    }

    if (modes.isLowSpeed) {
      signVc = x0_data[0];
    }
  }

  Vs_tmp = Vsy * (postProInputs->kappa + Vc);
  Vsy *= Vc;
  Vs_tmp = cos(tirParams->RCY1 * atan(Vs_tmp - Vsx * (Vs_tmp - atan(Vs_tmp)))) /
    cos(tirParams->RCY1 * atan(Vsy - Vsx * (Vsy - atan(Vsy)))) * Fy0 + signVc;
  Vsy = postProInputs->Fz / tirParams->FZMIN;
  signVc = postProInputs->Fz;
  if (postProInputs->Fz < tirParams->FZMIN) {
    signVc = postProInputs->Fz * (Vsy * Vsy);
  }

  Exa = fabs(postProInputs->gamma);
  Vsy = tirParams->QSX6 * signVc / tirParams->FNOMIN;
  forcesAndmoments->Mx = tirParams->UNLOADED_RADIUS * signVc * tirParams->LMX *
    (((((tirParams->QSX1 * tirParams->LVMX - tirParams->QSX2 *
         postProInputs->gamma * (tirParams->PPMX1 * incrVar_dpi + 1.0)) -
        tirParams->QSX12 * postProInputs->gamma * Exa) + tirParams->QSX3 *
       Vs_tmp / tirParams->FNOMIN) + tirParams->QSX4 * cos(tirParams->QSX5 *
       atan(Vsy * Vsy)) * sin(tirParams->QSX7 * postProInputs->gamma +
       tirParams->QSX8 * atan(tirParams->QSX9 * Vs_tmp / tirParams->FNOMIN))) +
     tirParams->QSX10 * atan(tirParams->QSX11 * signVc / tirParams->FNOMIN) *
     postProInputs->gamma) + tirParams->UNLOADED_RADIUS * Vs_tmp *
    tirParams->LMX * (tirParams->QSX13 + tirParams->QSX14 * Exa);
  signVc = postProInputs->uFz;
  if (postProInputs->uFz < tirParams->FZMIN) {
    signVc = postProInputs->uFz * (postProInputs->uFz / tirParams->FZMIN);
  }

  if (tirParams->FITTYP == 6.0) {
    Exa = postProInputs->uVcx / tirParams->LONGVL;
    Vc = -tirParams->UNLOADED_RADIUS * signVc * tirParams->LMY *
      (((tirParams->QSY1 + tirParams->QSY2 * (Fx / tirParams->FNOMIN)) +
        tirParams->QSY3 * fabs(Exa)) + tirParams->QSY4 * rt_powd_snf(Exa, 4.0));
  } else {
    Exa = postProInputs->uVcx / tirParams->LONGVL;
    Vsy = signVc / tirParams->FNOMIN;
    Vc = -tirParams->UNLOADED_RADIUS * tirParams->FNOMIN * tirParams->LMY *
      ((((tirParams->QSY1 + tirParams->QSY2 * (Fx / tirParams->FNOMIN)) +
         tirParams->QSY3 * fabs(Exa)) + tirParams->QSY4 * rt_powd_snf(Exa, 4.0))
       + (tirParams->QSY5 + tirParams->QSY6 * Vsy) * (postProInputs->gamma *
        postProInputs->gamma)) * (rt_powd_snf(Vsy, tirParams->QSY7) *
      rt_powd_snf(postProInputs->p / tirParams->NOMPRES, tirParams->QSY8));
  }

  loop_ub = 0;
  if (postProInputs->uVcx < 0.0) {
    loop_ub = 1;
  }

  if (0 <= loop_ub - 1) {
    x0_data[0] = -Vc;
  }

  signVc = Vc;
  if (postProInputs->uVcx < 0.0) {
    signVc = x0_data[0];
  }

  Vc = signVc;
  Vsy = tirParams->VXLOW / Vsx_tmp_tmp - 1.0;
  Exa = (-1.0 - tirParams->VXLOW) - Vsy;
  if ((postProInputs->ukappa >= Exa) && (postProInputs->ukappa <= Vsy)) {
    x0_data[0] = -1.0;

    /* INTERPOLATE Interpolate between two points. */
    /*  doExtras */
    for (i = 0; i < 1; i++) {
      x0_data[0] = (0.0 * (Vsy - postProInputs->ukappa) + 1.5707963267948966 *
                    (postProInputs->ukappa - x0_data[0])) / (Vsy - x0_data[0]);
    }

    for (loop_ub = 0; loop_ub < 1; loop_ub++) {
      x0_data[0] = sin(x0_data[0]);
    }

    for (i = 0; i < 1; i++) {
      x0_data[0] *= signVc;
    }

    signVc = x0_data[0];
    Vc = x0_data[0];
  }

  loop_ub = 0;
  if (postProInputs->ukappa < Exa) {
    loop_ub = 1;
  }

  if (0 <= loop_ub - 1) {
    x0_data[0] = -Vc;
  }

  if (postProInputs->ukappa < Exa) {
    signVc = x0_data[0];
  }

  if (postProInputs->uVcx == 0.0) {
    signVc = 0.0;
  }

  forcesAndmoments->My = signVc;
  j_expl_temp.LMUY_star = LMUY_star;
  j_expl_temp.LMUX_star = LMUX_star;
  j_expl_temp.gamma_star = starVar_gamma_star;
  j_expl_temp.alpha_star = starVar_alpha_star;
  k_expl_temp.LMUY_prime = primeVar_LMUY_prime;
  k_expl_temp.LMUX_prime = primeVar_LMUX_prime;
  k_expl_temp.alpha_prime = primeVar_alpha_prime;
  k_expl_temp.Fz0_prime = Fz0_prime;
  l_expl_temp.dpi = incrVar_dpi;
  l_expl_temp.dfz = incrVar_dfz;
  Solver_calculateMz(tirParams, postProInputs, &b_internalParams, modes,
                     j_expl_temp, k_expl_temp, l_expl_temp, alphar, alphat, Kxk,
                     Kya_prime, Vs_tmp, Fx, Dr, Cr, Br, Dt, Ct, Bt, Et,
                     &forcesAndmoments->Mz, &varinf->t, &varinf->Mzr);

  /*  Pack outputs */
  forcesAndmoments->Fx = Fx;
  forcesAndmoments->Fy = Vs_tmp;
  forcesAndmoments->Fz = postProInputs->Fz;
  varinf->Kxk = Kxk;
  varinf->Kya = Kya;
  varinf->muy = muy;
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const double inputs[6]
 *                double useMode
 *                struct_T *postProInputs
 *                h_struct_T *internalParams
 *                b_struct_T *modes
 * Return Type  : void
 */
void Solver_parseInputs(const struct0_T *tirParams, const double inputs[6],
  double useMode, struct_T *postProInputs, h_struct_T *internalParams,
  b_struct_T *modes)
{
  double hundreds;
  double tens;
  boolean_T useLimitsCheck;
  boolean_T useAlphaStar;
  boolean_T useTurnSlip;
  double alpha;
  double b_gamma;
  double phit;
  double dpi;
  double Fz;
  boolean_T isLowSpeed;
  int trueCount;
  int reductionSharp_size_idx_1;
  double reductionSharp_data[1];
  double reductionLinear_data[1];
  boolean_T isLowSpeedAlpha;
  int k;
  double reductionLinear_alpha_data[1];
  int i;
  double tmp_data[1];
  double unnamed_idx_0;
  double omega;
  hundreds = floor(useMode / 100.0);
  tens = floor((useMode - hundreds * 100.0) / 10.0);
  switch ((int)hundreds) {
   case 1:
    useLimitsCheck = true;
    break;

   case 2:
    useLimitsCheck = false;
    break;
  }

  switch ((int)tens) {
   case 1:
    useAlphaStar = false;
    break;

   case 2:
    useAlphaStar = true;
    break;
  }

  switch ((int)floor((useMode - hundreds * 100.0) - tens * 10.0)) {
   case 1:
    useTurnSlip = false;
    break;

   case 2:
    useTurnSlip = true;
    break;
  }

  postProInputs->kappa = inputs[1];
  alpha = inputs[2];
  b_gamma = inputs[3];
  phit = inputs[4];
  dpi = inputs[0];
  if (inputs[0] < 0.0) {
    dpi = 0.0;
  }

  Fz = dpi;
  postProInputs->uFz = dpi;
  postProInputs->ukappa = inputs[1];
  postProInputs->ukappaLow = inputs[1];
  tens = inputs[2];
  postProInputs->ugamma = inputs[3];
  postProInputs->uphit = inputs[4];
  hundreds = tirParams->INFLPRES;
  postProInputs->Fz_lowLimit = dpi;
  if (useLimitsCheck) {
    phit = inputs[4] * cos(inputs[2]);
    hundreds = fabs(inputs[5]);
    isLowSpeed = (hundreds <= tirParams->VXLOW);
    trueCount = 0;
    if (isLowSpeed) {
      trueCount = 1;
    }

    if (0 <= trueCount - 1) {
      reductionSharp_data[0] = 3.1415926535897931 * (inputs[5] /
        tirParams->VXLOW);
    }

    for (k = 0; k < trueCount; k++) {
      reductionSharp_data[0] = cos(reductionSharp_data[0]);
    }

    internalParams->reductionSmooth.size[0] = 1;
    internalParams->reductionSmooth.size[1] = trueCount;
    for (i = 0; i < trueCount; i++) {
      internalParams->reductionSmooth.data[0] = 1.0 - 0.5 *
        (reductionSharp_data[0] + 1.0);
    }

    trueCount = 0;
    if (isLowSpeed) {
      trueCount = 1;
    }

    if (0 <= trueCount - 1) {
      reductionSharp_data[0] = inputs[5] / tirParams->VXLOW;
    }

    if (0 <= trueCount - 1) {
      reductionLinear_data[0] = fabs(reductionSharp_data[0]);
    }

    reductionSharp_size_idx_1 = trueCount;
    if (0 <= trueCount - 1) {
      reductionSharp_data[0] = reductionLinear_data[0] * 1.5707963267948966;
    }

    for (k = 0; k < trueCount; k++) {
      reductionSharp_data[0] = sin(reductionSharp_data[0]);
    }

    if (0 <= trueCount - 1) {
      tmp_data[0] = inputs[1] * reductionLinear_data[0];
    }

    unnamed_idx_0 = inputs[1];
    if (isLowSpeed) {
      unnamed_idx_0 = tmp_data[0];
    }

    postProInputs->ukappaLow = unnamed_idx_0;
    if (0 <= trueCount - 1) {
      tmp_data[0] = phit * reductionLinear_data[0];
    }

    unnamed_idx_0 = phit;
    k = 0;
    if (isLowSpeed) {
      unnamed_idx_0 = tmp_data[0];
    }

    if (inputs[5] < 0.0) {
      k = 1;
    }

    if (0 <= k - 1) {
      tmp_data[0] = -unnamed_idx_0;
    }

    hundreds += fabs(tan(inputs[2]) * inputs[5]);
    isLowSpeedAlpha = (hundreds < tirParams->VXLOW);
    k = 0;
    if (inputs[5] < 0.0) {
      unnamed_idx_0 = tmp_data[0];
    }

    if (isLowSpeedAlpha) {
      k = 1;
    }

    phit = unnamed_idx_0;
    if (0 <= k - 1) {
      reductionLinear_alpha_data[0] = hundreds / tirParams->VXLOW;
    }

    if (0 <= trueCount - 1) {
      tmp_data[0] = inputs[1] * reductionLinear_data[0];
    }

    unnamed_idx_0 = inputs[1];
    if (isLowSpeed) {
      unnamed_idx_0 = tmp_data[0];
    }

    if (0 <= k - 1) {
      tmp_data[0] = inputs[2] * reductionLinear_alpha_data[0];
    }

    hundreds = inputs[2];
    if (isLowSpeedAlpha) {
      hundreds = tmp_data[0];
    }

    if (hundreds < tirParams->ALPMIN) {
      hundreds = tirParams->ALPMIN;
    }

    if (hundreds > tirParams->ALPMAX) {
      hundreds = tirParams->ALPMAX;
    }

    alpha = hundreds;
    hundreds = inputs[3];
    if (inputs[3] < tirParams->CAMMIN) {
      hundreds = tirParams->CAMMIN;
    }

    if (hundreds > tirParams->CAMMAX) {
      hundreds = tirParams->CAMMAX;
    }

    b_gamma = hundreds;
    hundreds = dpi;
    if (dpi > tirParams->FZMAX) {
      hundreds = tirParams->FZMAX;
    }

    Fz = hundreds;
    dpi = hundreds;
    if (hundreds < tirParams->FZMIN) {
      dpi = tirParams->FZMIN;
    }

    postProInputs->Fz_lowLimit = dpi;
    dpi = tirParams->INFLPRES;
    if (tirParams->INFLPRES < tirParams->PRESMIN) {
      dpi = tirParams->PRESMIN;
    }

    hundreds = dpi;
    if (dpi > tirParams->PRESMAX) {
      hundreds = tirParams->PRESMAX;
    }

    dpi = unnamed_idx_0;
    if (unnamed_idx_0 < tirParams->KPUMIN) {
      dpi = tirParams->KPUMIN;
    }

    unnamed_idx_0 = dpi;
    if (dpi > tirParams->KPUMAX) {
      unnamed_idx_0 = tirParams->KPUMAX;
    }

    postProInputs->kappa = unnamed_idx_0;
  } else {
    isLowSpeed = false;
    internalParams->reductionSmooth.size[0] = 1;
    internalParams->reductionSmooth.size[1] = 1;
    internalParams->reductionSmooth.data[0] = 1.0;
    reductionSharp_size_idx_1 = 1;
    reductionSharp_data[0] = 1.0;
    trueCount = 1;
    reductionLinear_data[0] = 1.0;
    isLowSpeedAlpha = false;
    k = 1;
    reductionLinear_alpha_data[0] = 1.0;
  }

  modes->useLimitsCheck = useLimitsCheck;
  modes->useAlphaStar = useAlphaStar;
  modes->useTurnSlip = useTurnSlip;
  modes->isLowSpeed = isLowSpeed;
  modes->isLowSpeedAlpha = isLowSpeedAlpha;
  modes->userDynamics = 0.0;
  internalParams->epsilonx = 1.0E-6;
  internalParams->epsilonk = 1.0E-6;
  internalParams->epsilony = 1.0E-6;
  internalParams->epsilonr = 1.0E-6;
  internalParams->epsilonv = 1.0E-6;
  internalParams->reductionSharp.size[0] = 1;
  internalParams->reductionSharp.size[1] = reductionSharp_size_idx_1;
  for (i = 0; i < reductionSharp_size_idx_1; i++) {
    internalParams->reductionSharp.data[0] = reductionSharp_data[0];
  }

  internalParams->reductionLinear.size[0] = 1;
  internalParams->reductionLinear.size[1] = trueCount;
  for (i = 0; i < trueCount; i++) {
    internalParams->reductionLinear.data[0] = reductionLinear_data[0];
  }

  internalParams->reductionLinear_alpha.size[0] = 1;
  internalParams->reductionLinear_alpha.size[1] = k;
  for (i = 0; i < k; i++) {
    internalParams->reductionLinear_alpha.data[0] = reductionLinear_alpha_data[0];
  }

  internalParams->zeta2 = 1.0;
  internalParams->epsilong = 0.0;
  postProInputs->alpha = alpha;
  postProInputs->gamma = b_gamma;
  postProInputs->phit = phit;
  postProInputs->Fz = Fz;
  postProInputs->p = hundreds;
  postProInputs->phi = 0.0;
  if (inputs[5] == 0.0) {
    tens = 0.0;
  }

  postProInputs->ualpha = tens;
  postProInputs->uVcx = inputs[5];
  postProInputs->nInputs = 6.0;
  dpi = (hundreds - tirParams->NOMPRES) / tirParams->NOMPRES;
  Solver_calculateRe(tirParams->LONGVL, tirParams->UNLOADED_RADIUS,
                     tirParams->FNOMIN, tirParams->VERTICAL_STIFFNESS,
                     tirParams->BREFF, tirParams->DREFF, tirParams->FREFF,
                     tirParams->Q_RE0, tirParams->Q_V1, tirParams->PFZ1, 0.0,
                     postProInputs->uFz, inputs[1], inputs[5], dpi, &hundreds,
                     &tens, &omega);
  postProInputs->omega = omega;
  Solver_calculateRe(tirParams->LONGVL, tirParams->UNLOADED_RADIUS,
                     tirParams->FNOMIN, tirParams->VERTICAL_STIFFNESS,
                     tirParams->BREFF, tirParams->DREFF, tirParams->FREFF,
                     tirParams->Q_RE0, tirParams->Q_V1, tirParams->PFZ1, omega,
                     postProInputs->uFz, inputs[1], inputs[5], dpi, &hundreds,
                     &tens, &unnamed_idx_0);
  postProInputs->Vsx = inputs[5] - hundreds * omega;
  if (useTurnSlip) {
    tens = tirParams->LFZO * tirParams->FNOMIN;
    unnamed_idx_0 = fabs(inputs[5]);
    hundreds = -tan(alpha) * unnamed_idx_0;
    dpi = sqrt(inputs[5] * inputs[5] + hundreds * hundreds);
    tens = tirParams->PECP1 * (tirParams->PECP2 * ((Fz - tens) / tens) + 1.0);
    isLowSpeed = (unnamed_idx_0 < tirParams->VXLOW);
    trueCount = -1;
    if (isLowSpeed) {
      trueCount = 0;
    }

    k = trueCount + 1;
    if (0 <= k - 1) {
      reductionSharp_data[0] = inputs[5];
    }

    for (k = 0; k <= trueCount; k++) {
      hundreds = reductionSharp_data[0];
      if (reductionSharp_data[0] < 0.0) {
        hundreds = -1.0;
      } else if (reductionSharp_data[0] > 0.0) {
        hundreds = 1.0;
      } else {
        if (reductionSharp_data[0] == 0.0) {
          hundreds = 0.0;
        }
      }

      reductionSharp_data[0] = hundreds;
    }

    for (k = 0; k <= trueCount; k++) {
      if (reductionSharp_data[0] == 0.0) {
        reductionSharp_data[0] = 1.0;
      }
    }

    for (i = 0; i <= trueCount; i++) {
      reductionSharp_data[0] *= tirParams->VXLOW;
    }

    if (isLowSpeed) {
      dpi = reductionSharp_data[0];
    }

    postProInputs->phi = 1.0 / dpi * (-phit * dpi - (1.0 - tens) * omega * sin
      (b_gamma));
    internalParams->epsilong = tens;
  }

  /*  costFz */
  /*  Protected methods */
  /*  Classdef */
}

/*
 * Arguments    : const struct0_T *tirParams
 *                const struct_T *postProInputs
 *                double omega
 *                double Romega
 *                double dpi
 *                const f_struct_T forcesAndmoments
 *                double Rl
 * Return Type  : double
 */
double __anon_fcn(const struct0_T *tirParams, const struct_T *postProInputs,
                  double omega, double Romega, double dpi, const f_struct_T
                  forcesAndmoments, double Rl)
{
  double a;
  double a_tmp;
  double rho_zfr;
  double b_a;
  double c_a;
  double rho_zg;
  double dv[1];
  int trueCount;

  /*  Cost function to calculate Rl. This measures the error */
  /*  between the tatget Fz (input of MFeval) vs the Fz calculated */
  /*  with the current Rl value. */
  /*  interpolate */
  /*  Calculate Fz from MF 6.2  with current Rl */
  /*  Cost (error) */
  a = tirParams->VERTICAL_STIFFNESS * tirParams->UNLOADED_RADIUS /
    tirParams->FNOMIN;
  a_tmp = Rl / Romega;
  rho_zfr = Romega - Rl;
  if (!(rho_zfr > 0.0)) {
    rho_zfr = 0.0;
  }

  b_a = (tirParams->Q_CAM1 * Rl + tirParams->Q_CAM2 * (Rl * Rl)) *
    postProInputs->gamma;
  c_a = (tirParams->Q_CAM1 * Romega + tirParams->Q_CAM2 * (Romega * Romega)) *
    postProInputs->gamma;
  rho_zg = b_a * b_a * ((1.075 - 0.5 * tirParams->ASPECT_RATIO) *
                        tirParams->WIDTH / 8.0) * fabs(tan(postProInputs->gamma))
    / (c_a * c_a) - tirParams->Q_CAM3 * rho_zfr * fabs(postProInputs->gamma);
  dv[0] = rho_zg;
  trueCount = 0;
  if (rtIsNaN(rho_zg)) {
    trueCount = 1;
  }

  if (0 <= trueCount - 1) {
    dv[0] = 0.0;
  }

  b_a = tirParams->Q_FCX * forcesAndmoments.Fx / tirParams->FNOMIN;
  rho_zg = rho_zfr + dv[0];
  if (!(rho_zg > 0.0)) {
    rho_zg = 0.0;
  }

  rho_zg /= tirParams->UNLOADED_RADIUS;
  c_a = rt_powd_snf(rho_zg, tirParams->Q_FCY2) * (tirParams->Q_FCY *
    (forcesAndmoments.Fy - ((tirParams->Q_FYS1 + tirParams->Q_FYS2 * a_tmp) +
    tirParams->Q_FYS3 * (a_tmp * a_tmp)) * postProInputs->gamma) /
    tirParams->FNOMIN);
  a = (((tirParams->Q_V2 * (tirParams->UNLOADED_RADIUS / tirParams->LONGVL) *
         fabs(omega) + 1.0) - b_a * b_a) - c_a * c_a) * (tirParams->PFZ1 * dpi +
    1.0) * (sqrt(a * a - 4.0 * tirParams->Q_FZ2) * rho_zg + tirParams->Q_FZ2 *
            (rho_zg * rho_zg)) * tirParams->FNOMIN - forcesAndmoments.Fz;
  return a * a;
}

/*
 * File trailer for Solver.c
 *
 * [EOF]
 */
