/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mfeval_wrapper.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

/* Include Files */
#include "mfeval_wrapper.h"
#include "Solver.h"
#include "fminsearch.h"
#include "mfeval_wrapper_data.h"
#include "mfeval_wrapper_initialize.h"
#include "mfeval_wrapper_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */

/*
 * Pack inputs
 * Arguments    : const struct0_T *parameterSource
 *                double Fz
 *                double kappa
 *                double alpha
 *                double b_gamma
 *                double phit
 *                double Vx
 *                double useMode
 *                double outMF[30]
 * Return Type  : void
 */
void mfeval_wrapper(const struct0_T *parameterSource, double Fz, double kappa,
                    double alpha, double b_gamma, double phit, double Vx, double
                    useMode, double outMF[30])
{
  double b_Fz[6];
  struct_T postProInputs;
  h_struct_T internalParams;
  b_struct_T modes;
  f_struct_T forcesAndmoments;
  g_struct_T varinf;
  double varinf_muy;
  double dpi_tmp;
  double Vcx;
  double Fz0;
  double R0;
  double omega;
  double a;
  double Romega;
  double Re;
  double Cz_tmp;
  double b_Cz_tmp;
  double Re_old;
  cell_0 fun_tunableEnvironment;
  double rho_zg;
  double dv[1];
  int trueCount;
  if (!isInitialized_mfeval_wrapper) {
    mfeval_wrapper_initialize();
  }

  /*  Call mfeval */
  /* MFEVAL evaluates Magic Formula 5.2, 6.1 or 6.2 tyre models in steady state */
  /* for series of input variables. */
  /*  */
  /*  outMF = mfeval(parameterSource, inputsMF, useMode) */
  /*  */
  /*  The formulation includes combined force/moment and turn slip */
  /*  calculations. */
  /*  ISO-W (TYDEX W) Contact-Patch Axis System coordinate system is used in */
  /*  all calculations. */
  /*  All the units will be SI (N,m,s,rad,kg) */
  /*  */
  /*  parameterSource refers to a MF-Tyre tyre property file (.TIR) containing */
  /*  all the Magic Formula coefficients or to a structure with all the */
  /*  parameters. */
  /*  */
  /*  inputsMF = [Fz kappa alpha gamma phit Vx P* omega*], where */
  /*  Fz     = normal load on the tyre  [N] */
  /*  kappa  = longitudinal slip        [dimensionless, -1: locked wheel] */
  /*  alpha  = side slip angle          [rad] */
  /*  gamma  = inclination angle        [rad] */
  /*  phit   = turn slip                [1/m] */
  /*  Vx     = forward velocity         [m/s] */
  /*  P*     = pressure                 [Pa] */
  /*  omega* = rotational speed         [rad/s] */
  /*  */
  /*  P* and omega* are optional inputs. If they are not included pressure is */
  /*  constant and equal to the inflation pressure on the TIR file and the */
  /*  rotational speed is approximated. */
  /*  */
  /*  useMode specifies the type of calculation performed: */
  /*       1: combined force/moment calculation */
  /*       2: combined force/moment calculation + turn slip */
  /*     +10: revoke alpha_star definition */
  /*     +20: include alpha_star definition */
  /*    +100: include limit checks */
  /*    +200: ignore limit checks */
  /*  */
  /*  For example: useMode = 121 implies: */
  /*    -combined force/moment */
  /*    -include alpha_star */
  /*    -include limit checks */
  /*  */
  /*  For normal situations turn slip may be neglected, meaning that the radius */
  /*  of the path is close to infinity (straight line). */
  /*  Alpha_star improves the accuracy of the model for very large slip angles */
  /*  and possibly backward running of the wheel. */
  /*  The limit checks verify that the inputs are inside the stable range of */
  /*  the model. */
  /*  */
  /*  outMF consists of 30 columns: */
  /*  1 - Fx: longitudinal force         16 - t: pneumatic trail */
  /*  2 - Fy: lateral force              17 - mux: longitudinal friction coeff. */
  /*  3 - Fz: normal force               18 - muy: lateral friction coeff. */
  /*  4 - Mx: overturning moment         19 - omega: rot. speed */
  /*  5 - My: rolling resistance moment  20 - Rl : loaded radius */
  /*  6 - Mz: self aligning moment       21 - 2b: contact patch width */
  /*  7 - kappa: longitudinal slip       22 - Mzr: residual torque */
  /*  8 - alpha: side slip angle         23 - Cx: longitudinal stiffness */
  /*  9 - gamma: inclination angle       24 - Cy: lateral stiffness */
  /*  10 - phit: turn slip               25 - Cz: vertical stiffness */
  /*  11 - Vx: longitudinal velocity     26 - Kya: cornering stiffness */
  /*  12 - P: pressure                   27 - sigmax: long. relax. length */
  /*  13 - Re: effective rolling radius  28 - sigmay: lat. relax. length */
  /*  14 - rho: tyre deflection          29 - iKya: Instantaneous cornering stiff. */
  /*  15 - 2a: contact patch length      30 - Kxk: slip stiffness */
  /*  */
  /*  The equations here presented are published in the book: */
  /*  Title:    	Tire and Vehicle Dynamics */
  /*  Author:       Hans Pacejka */
  /*  Edition:      3, revised */
  /*  Publisher:	Elsevier, 2012 */
  /*  ISBN:         0080970176, 9780080970172 */
  /*  Length:       672 pages */
  /*  */
  /*  And in the following paper: */
  /*  Besselink, I. J. M. , Schmeitz, A. J. C. and Pacejka, H. B.(2010) 'An */
  /*  improved Magic Formula/Swift tyre model that can handle inflation */
  /*  pressure changes', Vehicle System Dynamics, 48: 1, 337 — 352 */
  /*  DOI: 10.1080/00423111003748088 */
  /*  URL: http://dx.doi.org/10.1080/00423111003748088 */
  /*  */
  /*  Code developed by:	Marco Furlan */
  /*  Email address:        marcofurlan92@gmail.com */
  /*  */
  /*  See the <a href="matlab:web(fullfile(mfevalroot, '..', 'doc','index.html'))">documentation</a> for more details and examples */
  /*  Versions change log: */
  /*  Please refer to the documentation: */
  /*  "help mfeval" > click on the documentation hyperlink > Release Notes */
  /*  Declare extrinsic functions for C code generation */
  /*  Check number of inputs and assign */
  /*  Validate parameterSource */
  /*  Validate inputsMF */
  /*  Check the number of columns of inputsMF */
  /*  Validate useMode */
  /*  Check the number of digits of useMode */
  /*  Load the TIR file parameters */
  /*  Structure of parameter */
  /*  If else */
  /*  Validate the Magic Formula Version */
  /*  Run solver */
  b_Fz[0] = Fz;
  b_Fz[1] = kappa;
  b_Fz[2] = alpha;
  b_Fz[3] = b_gamma;
  b_Fz[4] = phit;
  b_Fz[5] = Vx;
  Solver_parseInputs(parameterSource, b_Fz, useMode, &postProInputs,
                     &internalParams, &modes);
  Solver_doForcesAndMoments(parameterSource, &postProInputs, &internalParams,
    modes, &forcesAndmoments, &varinf);
  varinf_muy = varinf.muy;

  /*  doForcesAndMoments */
  dpi_tmp = (postProInputs.p - parameterSource->NOMPRES) /
    parameterSource->NOMPRES;

  /*  [Eqn (4.E2b) Page 177 - Book] */
  Vcx = postProInputs.uVcx;
  Fz0 = parameterSource->FNOMIN;
  R0 = parameterSource->UNLOADED_RADIUS;
  omega = postProInputs.omega;
  a = postProInputs.omega * parameterSource->UNLOADED_RADIUS /
    parameterSource->LONGVL;
  Romega = parameterSource->UNLOADED_RADIUS * (parameterSource->Q_RE0 +
    parameterSource->Q_V1 * (a * a));
  Re = parameterSource->UNLOADED_RADIUS * 0.965;
  Cz_tmp = parameterSource->PFZ1 * dpi_tmp + 1.0;
  b_Cz_tmp = parameterSource->VERTICAL_STIFFNESS * Cz_tmp;
  Re_old = parameterSource->UNLOADED_RADIUS;
  while (fabs(Re_old - Re) > 1.0E-9) {
    Re_old = Re;
    omega = (postProInputs.ukappa * Vcx + Vcx) / Re;
    a = omega * R0 / parameterSource->LONGVL;
    Romega = R0 * (parameterSource->Q_RE0 + parameterSource->Q_V1 * (a * a));
    rho_zg = postProInputs.uFz / Fz0;
    Re = Romega - Fz0 / b_Cz_tmp * (parameterSource->DREFF * atan
      (parameterSource->BREFF * rho_zg) + parameterSource->FREFF * rho_zg);
  }

  /*  Calculate the radius with the MF61 or MF62 formulation */
  if (parameterSource->FITTYP == 61.0) {
    a = parameterSource->VERTICAL_STIFFNESS * parameterSource->UNLOADED_RADIUS /
      parameterSource->FNOMIN;
    rho_zg = 4.0 * parameterSource->Q_FZ2;
    Fz0 = sqrt(a * a - rho_zg);
    a = parameterSource->Q_FCX * forcesAndmoments.Fx / parameterSource->FNOMIN;
    Vcx = parameterSource->Q_FCY * forcesAndmoments.Fy / parameterSource->FNOMIN;
    Vcx = -(postProInputs.Fz_lowLimit / ((((parameterSource->Q_V2 *
                (parameterSource->UNLOADED_RADIUS / parameterSource->LONGVL) *
                fabs(omega) + 1.0) - a * a) - Vcx * Vcx) * Cz_tmp *
             parameterSource->FNOMIN));
    rho_zg = Fz0 * Fz0 - rho_zg * Vcx;
    if (rho_zg > 0.0) {
      rho_zg = (-Fz0 + sqrt(Fz0 * Fz0 - 4.0 * parameterSource->Q_FZ2 * Vcx)) /
        (2.0 * parameterSource->Q_FZ2);
    } else {
      rho_zg = (-Fz0 + sqrt(rho_zg)) / (2.0 * parameterSource->Q_FZ2);
    }

    Fz0 = rho_zg * parameterSource->UNLOADED_RADIUS;
    if (!(Fz0 > 0.0)) {
      Fz0 = 0.0;
    }

    Cz_tmp = Romega - Fz0;
    if (!(Cz_tmp > 0.0)) {
      Cz_tmp = 0.0;
    }

    dv[0] = Fz0;
    if (postProInputs.Fz_lowLimit == 0.0) {
      dv[0] = 1.0E-6;
    }

    outMF[24] = postProInputs.Fz_lowLimit / dv[0];
    outMF[13] = dv[0];
  } else {
    fun_tunableEnvironment.f4 = omega;
    fun_tunableEnvironment.f5 = Romega;
    fun_tunableEnvironment.f6 = dpi_tmp;
    fun_tunableEnvironment.f2 = *parameterSource;
    fun_tunableEnvironment.f3 = postProInputs;
    fun_tunableEnvironment.f7 = forcesAndmoments;
    Cz_tmp = fminsearch(&fun_tunableEnvironment);
    Fz0 = Romega - Cz_tmp;
    if (!(Fz0 > 0.0)) {
      Fz0 = 0.0;
    }

    a = (parameterSource->Q_CAM1 * Cz_tmp + parameterSource->Q_CAM2 * (Cz_tmp *
          Cz_tmp)) * postProInputs.gamma;
    Vcx = (parameterSource->Q_CAM1 * Romega + parameterSource->Q_CAM2 * (Romega *
            Romega)) * postProInputs.gamma;
    rho_zg = a * a * ((1.075 - 0.5 * parameterSource->ASPECT_RATIO) *
                      parameterSource->WIDTH / 8.0) * fabs(tan
      (postProInputs.gamma)) / (Vcx * Vcx) - parameterSource->Q_CAM3 * Fz0 *
      fabs(postProInputs.gamma);
    dv[0] = rho_zg;
    trueCount = 0;
    if (rtIsNaN(rho_zg)) {
      trueCount = 1;
    }

    if (0 <= trueCount - 1) {
      dv[0] = 0.0;
    }

    rho_zg = Fz0 + dv[0];
    if (rho_zg > 0.0) {
      dv[0] = rho_zg;
    } else {
      dv[0] = 0.0;
    }

    if (forcesAndmoments.Fz == 0.0) {
      dv[0] = 1.0E-6;
    }

    outMF[13] = dv[0];
    outMF[24] = forcesAndmoments.Fz / dv[0];
  }

  dv[0] = postProInputs.uFz;
  if ((parameterSource->UNLOADED_RADIUS - postProInputs.uFz / b_Cz_tmp) -
      (parameterSource->RIM_RADIUS + parameterSource->BOTTOM_OFFST) < 0.0) {
    dv[0] = ((parameterSource->UNLOADED_RADIUS - parameterSource->RIM_RADIUS) -
             parameterSource->BOTTOM_OFFST) * b_Cz_tmp;
  }

  if (parameterSource->FITTYP == 6.0) {
    rho_zg = parameterSource->Q_A1;
    Vcx = parameterSource->Q_A2;
    if ((parameterSource->Q_A1 == 0.0) && (parameterSource->Q_A2 == 0.0)) {
      rho_zg = log10(parameterSource->UNLOADED_RADIUS *
                     (parameterSource->VERTICAL_STIFFNESS /
                      parameterSource->FNOMIN));
      rho_zg = ((-0.0388 * rt_powd_snf(rho_zg, 3.0) + 0.2509 * (rho_zg * rho_zg))
                + -0.6283 * rho_zg) + 0.6279;
      Vcx = 1.693 * (rho_zg * rho_zg);
    }

    Fz0 = dv[0] / parameterSource->FNOMIN;
    a = parameterSource->UNLOADED_RADIUS * (Vcx * Fz0 + rho_zg * sqrt(Fz0));
    Vcx = parameterSource->WIDTH / 2.0;
  } else {
    Fz0 = dv[0] / (b_Cz_tmp * parameterSource->UNLOADED_RADIUS);
    a = parameterSource->UNLOADED_RADIUS * (parameterSource->Q_RA2 * Fz0 +
      parameterSource->Q_RA1 * sqrt(Fz0));
    Vcx = parameterSource->WIDTH * (parameterSource->Q_RB2 * Fz0 +
      parameterSource->Q_RB1 * rt_powd_snf(Fz0, 0.33333333333333331));
  }

  Fz0 = parameterSource->LFZO * parameterSource->FNOMIN;
  R0 = (postProInputs.Fz - Fz0) / Fz0;
  rho_zg = R0 * R0;
  Re_old = parameterSource->LONGITUDINAL_STIFFNESS * ((parameterSource->PCFX1 *
    R0 + 1.0) + parameterSource->PCFX2 * rho_zg) * (parameterSource->PCFX3 *
    dpi_tmp + 1.0);
  rho_zg = parameterSource->LATERAL_STIFFNESS * ((parameterSource->PCFY1 * R0 +
    1.0) + parameterSource->PCFY2 * rho_zg) * (parameterSource->PCFY3 * dpi_tmp
    + 1.0);
  if (parameterSource->FITTYP == 6.0) {
    outMF[26] = (parameterSource->PTX1 + parameterSource->PTX2 * R0) * exp
      (-parameterSource->PTX3 * R0) * parameterSource->LSGKP *
      parameterSource->UNLOADED_RADIUS * postProInputs.Fz /
      parameterSource->FNOMIN;
    outMF[27] = parameterSource->PTY1 * sin(2.0 * atan(postProInputs.Fz /
      (parameterSource->PTY2 * Fz0))) * (1.0 - parameterSource->PKY3 * fabs
      (postProInputs.gamma)) * parameterSource->UNLOADED_RADIUS *
      parameterSource->LFZO * parameterSource->LSGAL;
  } else {
    outMF[26] = fabs(varinf.Kxk / Re_old);
    outMF[27] = fabs(varinf.Kya / rho_zg);
  }

  outMF[12] = Re;
  outMF[18] = omega;
  outMF[19] = Cz_tmp;
  outMF[22] = Re_old;
  outMF[23] = rho_zg;
  if (modes.useLimitsCheck && (varinf.muy < 0.0)) {
    varinf_muy = fabs(varinf.muy);
  }

  outMF[0] = forcesAndmoments.Fx;
  outMF[1] = forcesAndmoments.Fy;
  outMF[2] = postProInputs.uFz;
  outMF[3] = forcesAndmoments.Mx;
  outMF[4] = forcesAndmoments.My;
  outMF[5] = forcesAndmoments.Mz;
  outMF[6] = postProInputs.ukappaLow;
  outMF[7] = postProInputs.ualpha;
  outMF[8] = postProInputs.ugamma;
  outMF[9] = postProInputs.phit;
  outMF[10] = postProInputs.uVcx;
  outMF[11] = postProInputs.p;
  outMF[14] = 2.0 * a;
  outMF[15] = varinf.t;
  outMF[16] = varinf.mux;
  outMF[17] = varinf_muy;
  outMF[20] = 2.0 * Vcx;
  outMF[21] = varinf.Mzr;
  outMF[25] = varinf.Kya;
  outMF[28] = 0.0;
  outMF[29] = varinf.Kxk;

  /*  mfeval */
}

/*
 * File trailer for mfeval_wrapper.c
 *
 * [EOF]
 */
