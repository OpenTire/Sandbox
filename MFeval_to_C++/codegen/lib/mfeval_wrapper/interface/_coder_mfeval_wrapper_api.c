/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_mfeval_wrapper_api.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:15:56
 */

/* Include Files */
#include "_coder_mfeval_wrapper_api.h"
#include "_coder_mfeval_wrapper_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131594U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "mfeval_wrapper",                    /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[5]);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[7]);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *parameterSource, const char_T *identifier, struct0_T *y);
static const mxArray *emlrt_marshallOut(const real_T u[30]);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[8]);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[9]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[4]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[6]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[39]);
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Fz, const
  char_T *identifier);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[5]);
static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[7]);
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[8]);
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[9]);
static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[4]);
static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[6]);
static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[39]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct0_T *y
 * Return Type  : void
 */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[284] = { "FILE_TYPE", "FILE_VERSION",
    "FILE_FORMAT", "LENGTH", "FORCE", "ANGLE", "MASS", "TIME", "FITTYP",
    "TYRESIDE", "LONGVL", "VXLOW", "ROAD_INCREMENT", "ROAD_DIRECTION",
    "PROPERTY_FILE_FORMAT", "USER_SUB_ID", "N_TIRE_STATES", "USE_MODE",
    "HMAX_LOCAL", "TIME_SWITCH_INTEG", "UNLOADED_RADIUS", "WIDTH",
    "ASPECT_RATIO", "RIM_RADIUS", "RIM_WIDTH", "INFLPRES", "NOMPRES", "MASS1",
    "IXX", "IYY", "BELT_MASS", "BELT_IXX", "BELT_IYY", "GRAVITY", "FNOMIN",
    "VERTICAL_STIFFNESS", "VERTICAL_DAMPING", "MC_CONTOUR_A", "MC_CONTOUR_B",
    "BREFF", "DREFF", "FREFF", "Q_RE0", "Q_V1", "Q_V2", "Q_FZ2", "Q_FCX",
    "Q_FCY", "Q_CAM", "PFZ1", "BOTTOM_OFFST", "BOTTOM_STIFF",
    "LONGITUDINAL_STIFFNESS", "LATERAL_STIFFNESS", "YAW_STIFFNESS", "FREQ_LONG",
    "FREQ_LAT", "FREQ_YAW", "FREQ_WINDUP", "DAMP_LONG", "DAMP_LAT", "DAMP_YAW",
    "DAMP_WINDUP", "DAMP_RESIDUAL", "DAMP_VLOW", "Q_BVX", "Q_BVT", "PCFX1",
    "PCFX2", "PCFX3", "PCFY1", "PCFY2", "PCFY3", "PCMZ1", "Q_RA1", "Q_RA2",
    "Q_RB1", "Q_RB2", "ELLIPS_SHIFT", "ELLIPS_LENGTH", "ELLIPS_HEIGHT",
    "ELLIPS_ORDER", "ELLIPS_MAX_STEP", "ELLIPS_NWIDTH", "ELLIPS_NLENGTH",
    "PRESMIN", "PRESMAX", "FZMIN", "FZMAX", "KPUMIN", "KPUMAX", "ALPMIN",
    "ALPMAX", "CAMMIN", "CAMMAX", "LFZO", "LCX", "LMUX", "LEX", "LKX", "LHX",
    "LVX", "LCY", "LMUY", "LEY", "LKY", "LHY", "LVY", "LTR", "LRES", "LXAL",
    "LYKA", "LVYKA", "LS", "LKYC", "LKZC", "LVMX", "LMX", "LMY", "LMP", "PCX1",
    "PDX1", "PDX2", "PDX3", "PEX1", "PEX2", "PEX3", "PEX4", "PKX1", "PKX2",
    "PKX3", "PHX1", "PHX2", "PVX1", "PVX2", "PPX1", "PPX2", "PPX3", "PPX4",
    "RBX1", "RBX2", "RBX3", "RCX1", "REX1", "REX2", "RHX1", "QSX1", "QSX2",
    "QSX3", "QSX4", "QSX5", "QSX6", "QSX7", "QSX8", "QSX9", "QSX10", "QSX11",
    "QSX12", "QSX13", "QSX14", "PPMX1", "PCY1", "PDY1", "PDY2", "PDY3", "PEY1",
    "PEY2", "PEY3", "PEY4", "PEY5", "PKY1", "PKY2", "PKY3", "PKY4", "PKY5",
    "PKY6", "PKY7", "PHY1", "PHY2", "PVY1", "PVY2", "PVY3", "PVY4", "PPY1",
    "PPY2", "PPY3", "PPY4", "PPY5", "RBY1", "RBY2", "RBY3", "RBY4", "RCY1",
    "REY1", "REY2", "RHY1", "RHY2", "RVY1", "RVY2", "RVY3", "RVY4", "RVY5",
    "RVY6", "QSY1", "QSY2", "QSY3", "QSY4", "QSY5", "QSY6", "QSY7", "QSY8",
    "QBZ1", "QBZ2", "QBZ3", "QBZ4", "QBZ5", "QBZ9", "QBZ10", "QCZ1", "QDZ1",
    "QDZ2", "QDZ3", "QDZ4", "QDZ6", "QDZ7", "QDZ8", "QDZ9", "QDZ10", "QDZ11",
    "QEZ1", "QEZ2", "QEZ3", "QEZ4", "QEZ5", "QHZ1", "QHZ2", "QHZ3", "QHZ4",
    "PPZ1", "PPZ2", "SSZ1", "SSZ2", "SSZ3", "SSZ4", "PDXP1", "PDXP2", "PDXP3",
    "PKYP1", "PDYP1", "PDYP2", "PDYP3", "PDYP4", "PHYP1", "PHYP2", "PHYP3",
    "PHYP4", "PECP1", "PECP2", "QDTP1", "QCRP1", "QCRP2", "QBRP1", "QDRP1",
    "FUNCTION_NAME", "SWITCH_INTEG", "Q_FCY2", "Q_CAM1", "Q_CAM2", "Q_CAM3",
    "Q_FYS1", "Q_FYS2", "Q_FYS3", "ENV_C1", "ENV_C2", "Q_A1", "Q_A2", "PHY3",
    "PTX1", "PTX2", "PTX3", "PTY1", "PTY2", "LSGKP", "LSGAL" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 284, fieldNames, 0U, &dims);
  thisId.fIdentifier = "FILE_TYPE";
  c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 0, "FILE_TYPE")),
                     &thisId, y->FILE_TYPE);
  thisId.fIdentifier = "FILE_VERSION";
  y->FILE_VERSION = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 1, "FILE_VERSION")), &thisId);
  thisId.fIdentifier = "FILE_FORMAT";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 2,
    "FILE_FORMAT")), &thisId, y->FILE_FORMAT);
  thisId.fIdentifier = "LENGTH";
  e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3, "LENGTH")),
                     &thisId, y->LENGTH);
  thisId.fIdentifier = "FORCE";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4, "FORCE")),
                     &thisId, y->FORCE);
  thisId.fIdentifier = "ANGLE";
  g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5, "ANGLE")),
                     &thisId, y->ANGLE);
  thisId.fIdentifier = "MASS";
  h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6, "MASS")),
                     &thisId, y->MASS);
  thisId.fIdentifier = "TIME";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7, "TIME")),
                     &thisId, y->TIME);
  thisId.fIdentifier = "FITTYP";
  y->FITTYP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "FITTYP")), &thisId);
  thisId.fIdentifier = "TYRESIDE";
  i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9, "TYRESIDE")),
                     &thisId, y->TYRESIDE);
  thisId.fIdentifier = "LONGVL";
  y->LONGVL = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "LONGVL")), &thisId);
  thisId.fIdentifier = "VXLOW";
  y->VXLOW = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11,
    "VXLOW")), &thisId);
  thisId.fIdentifier = "ROAD_INCREMENT";
  y->ROAD_INCREMENT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 12, "ROAD_INCREMENT")), &thisId);
  thisId.fIdentifier = "ROAD_DIRECTION";
  y->ROAD_DIRECTION = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 13, "ROAD_DIRECTION")), &thisId);
  thisId.fIdentifier = "PROPERTY_FILE_FORMAT";
  i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14,
    "PROPERTY_FILE_FORMAT")), &thisId, y->PROPERTY_FILE_FORMAT);
  thisId.fIdentifier = "USER_SUB_ID";
  y->USER_SUB_ID = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 15, "USER_SUB_ID")), &thisId);
  thisId.fIdentifier = "N_TIRE_STATES";
  y->N_TIRE_STATES = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 16, "N_TIRE_STATES")), &thisId);
  thisId.fIdentifier = "USE_MODE";
  y->USE_MODE = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    17, "USE_MODE")), &thisId);
  thisId.fIdentifier = "HMAX_LOCAL";
  y->HMAX_LOCAL = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    18, "HMAX_LOCAL")), &thisId);
  thisId.fIdentifier = "TIME_SWITCH_INTEG";
  y->TIME_SWITCH_INTEG = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 19, "TIME_SWITCH_INTEG")), &thisId);
  thisId.fIdentifier = "UNLOADED_RADIUS";
  y->UNLOADED_RADIUS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 20, "UNLOADED_RADIUS")), &thisId);
  thisId.fIdentifier = "WIDTH";
  y->WIDTH = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 21,
    "WIDTH")), &thisId);
  thisId.fIdentifier = "ASPECT_RATIO";
  y->ASPECT_RATIO = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 22, "ASPECT_RATIO")), &thisId);
  thisId.fIdentifier = "RIM_RADIUS";
  y->RIM_RADIUS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    23, "RIM_RADIUS")), &thisId);
  thisId.fIdentifier = "RIM_WIDTH";
  y->RIM_WIDTH = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    24, "RIM_WIDTH")), &thisId);
  thisId.fIdentifier = "INFLPRES";
  y->INFLPRES = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    25, "INFLPRES")), &thisId);
  thisId.fIdentifier = "NOMPRES";
  y->NOMPRES = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    26, "NOMPRES")), &thisId);
  thisId.fIdentifier = "MASS1";
  y->MASS1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 27,
    "MASS1")), &thisId);
  thisId.fIdentifier = "IXX";
  y->IXX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 28,
    "IXX")), &thisId);
  thisId.fIdentifier = "IYY";
  y->IYY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 29,
    "IYY")), &thisId);
  thisId.fIdentifier = "BELT_MASS";
  y->BELT_MASS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    30, "BELT_MASS")), &thisId);
  thisId.fIdentifier = "BELT_IXX";
  y->BELT_IXX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    31, "BELT_IXX")), &thisId);
  thisId.fIdentifier = "BELT_IYY";
  y->BELT_IYY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    32, "BELT_IYY")), &thisId);
  thisId.fIdentifier = "GRAVITY";
  y->GRAVITY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    33, "GRAVITY")), &thisId);
  thisId.fIdentifier = "FNOMIN";
  y->FNOMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 34,
    "FNOMIN")), &thisId);
  thisId.fIdentifier = "VERTICAL_STIFFNESS";
  y->VERTICAL_STIFFNESS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 35, "VERTICAL_STIFFNESS")), &thisId);
  thisId.fIdentifier = "VERTICAL_DAMPING";
  y->VERTICAL_DAMPING = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 36, "VERTICAL_DAMPING")), &thisId);
  thisId.fIdentifier = "MC_CONTOUR_A";
  y->MC_CONTOUR_A = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 37, "MC_CONTOUR_A")), &thisId);
  thisId.fIdentifier = "MC_CONTOUR_B";
  y->MC_CONTOUR_B = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 38, "MC_CONTOUR_B")), &thisId);
  thisId.fIdentifier = "BREFF";
  y->BREFF = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 39,
    "BREFF")), &thisId);
  thisId.fIdentifier = "DREFF";
  y->DREFF = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 40,
    "DREFF")), &thisId);
  thisId.fIdentifier = "FREFF";
  y->FREFF = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 41,
    "FREFF")), &thisId);
  thisId.fIdentifier = "Q_RE0";
  y->Q_RE0 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 42,
    "Q_RE0")), &thisId);
  thisId.fIdentifier = "Q_V1";
  y->Q_V1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 43,
    "Q_V1")), &thisId);
  thisId.fIdentifier = "Q_V2";
  y->Q_V2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 44,
    "Q_V2")), &thisId);
  thisId.fIdentifier = "Q_FZ2";
  y->Q_FZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 45,
    "Q_FZ2")), &thisId);
  thisId.fIdentifier = "Q_FCX";
  y->Q_FCX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 46,
    "Q_FCX")), &thisId);
  thisId.fIdentifier = "Q_FCY";
  y->Q_FCY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 47,
    "Q_FCY")), &thisId);
  thisId.fIdentifier = "Q_CAM";
  y->Q_CAM = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 48,
    "Q_CAM")), &thisId);
  thisId.fIdentifier = "PFZ1";
  y->PFZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 49,
    "PFZ1")), &thisId);
  thisId.fIdentifier = "BOTTOM_OFFST";
  y->BOTTOM_OFFST = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 50, "BOTTOM_OFFST")), &thisId);
  thisId.fIdentifier = "BOTTOM_STIFF";
  y->BOTTOM_STIFF = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 51, "BOTTOM_STIFF")), &thisId);
  thisId.fIdentifier = "LONGITUDINAL_STIFFNESS";
  y->LONGITUDINAL_STIFFNESS = d_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetFieldR2017b(sp, u, 0, 52, "LONGITUDINAL_STIFFNESS")), &thisId);
  thisId.fIdentifier = "LATERAL_STIFFNESS";
  y->LATERAL_STIFFNESS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 53, "LATERAL_STIFFNESS")), &thisId);
  thisId.fIdentifier = "YAW_STIFFNESS";
  y->YAW_STIFFNESS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 54, "YAW_STIFFNESS")), &thisId);
  thisId.fIdentifier = "FREQ_LONG";
  y->FREQ_LONG = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    55, "FREQ_LONG")), &thisId);
  thisId.fIdentifier = "FREQ_LAT";
  y->FREQ_LAT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    56, "FREQ_LAT")), &thisId);
  thisId.fIdentifier = "FREQ_YAW";
  y->FREQ_YAW = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    57, "FREQ_YAW")), &thisId);
  thisId.fIdentifier = "FREQ_WINDUP";
  y->FREQ_WINDUP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 58, "FREQ_WINDUP")), &thisId);
  thisId.fIdentifier = "DAMP_LONG";
  y->DAMP_LONG = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    59, "DAMP_LONG")), &thisId);
  thisId.fIdentifier = "DAMP_LAT";
  y->DAMP_LAT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    60, "DAMP_LAT")), &thisId);
  thisId.fIdentifier = "DAMP_YAW";
  y->DAMP_YAW = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    61, "DAMP_YAW")), &thisId);
  thisId.fIdentifier = "DAMP_WINDUP";
  y->DAMP_WINDUP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 62, "DAMP_WINDUP")), &thisId);
  thisId.fIdentifier = "DAMP_RESIDUAL";
  y->DAMP_RESIDUAL = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 63, "DAMP_RESIDUAL")), &thisId);
  thisId.fIdentifier = "DAMP_VLOW";
  y->DAMP_VLOW = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    64, "DAMP_VLOW")), &thisId);
  thisId.fIdentifier = "Q_BVX";
  y->Q_BVX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 65,
    "Q_BVX")), &thisId);
  thisId.fIdentifier = "Q_BVT";
  y->Q_BVT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 66,
    "Q_BVT")), &thisId);
  thisId.fIdentifier = "PCFX1";
  y->PCFX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 67,
    "PCFX1")), &thisId);
  thisId.fIdentifier = "PCFX2";
  y->PCFX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 68,
    "PCFX2")), &thisId);
  thisId.fIdentifier = "PCFX3";
  y->PCFX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 69,
    "PCFX3")), &thisId);
  thisId.fIdentifier = "PCFY1";
  y->PCFY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 70,
    "PCFY1")), &thisId);
  thisId.fIdentifier = "PCFY2";
  y->PCFY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 71,
    "PCFY2")), &thisId);
  thisId.fIdentifier = "PCFY3";
  y->PCFY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 72,
    "PCFY3")), &thisId);
  thisId.fIdentifier = "PCMZ1";
  y->PCMZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 73,
    "PCMZ1")), &thisId);
  thisId.fIdentifier = "Q_RA1";
  y->Q_RA1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 74,
    "Q_RA1")), &thisId);
  thisId.fIdentifier = "Q_RA2";
  y->Q_RA2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 75,
    "Q_RA2")), &thisId);
  thisId.fIdentifier = "Q_RB1";
  y->Q_RB1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 76,
    "Q_RB1")), &thisId);
  thisId.fIdentifier = "Q_RB2";
  y->Q_RB2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 77,
    "Q_RB2")), &thisId);
  thisId.fIdentifier = "ELLIPS_SHIFT";
  y->ELLIPS_SHIFT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 78, "ELLIPS_SHIFT")), &thisId);
  thisId.fIdentifier = "ELLIPS_LENGTH";
  y->ELLIPS_LENGTH = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 79, "ELLIPS_LENGTH")), &thisId);
  thisId.fIdentifier = "ELLIPS_HEIGHT";
  y->ELLIPS_HEIGHT = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 80, "ELLIPS_HEIGHT")), &thisId);
  thisId.fIdentifier = "ELLIPS_ORDER";
  y->ELLIPS_ORDER = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 81, "ELLIPS_ORDER")), &thisId);
  thisId.fIdentifier = "ELLIPS_MAX_STEP";
  y->ELLIPS_MAX_STEP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 82, "ELLIPS_MAX_STEP")), &thisId);
  thisId.fIdentifier = "ELLIPS_NWIDTH";
  y->ELLIPS_NWIDTH = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 83, "ELLIPS_NWIDTH")), &thisId);
  thisId.fIdentifier = "ELLIPS_NLENGTH";
  y->ELLIPS_NLENGTH = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 84, "ELLIPS_NLENGTH")), &thisId);
  thisId.fIdentifier = "PRESMIN";
  y->PRESMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    85, "PRESMIN")), &thisId);
  thisId.fIdentifier = "PRESMAX";
  y->PRESMAX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    86, "PRESMAX")), &thisId);
  thisId.fIdentifier = "FZMIN";
  y->FZMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 87,
    "FZMIN")), &thisId);
  thisId.fIdentifier = "FZMAX";
  y->FZMAX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 88,
    "FZMAX")), &thisId);
  thisId.fIdentifier = "KPUMIN";
  y->KPUMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 89,
    "KPUMIN")), &thisId);
  thisId.fIdentifier = "KPUMAX";
  y->KPUMAX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 90,
    "KPUMAX")), &thisId);
  thisId.fIdentifier = "ALPMIN";
  y->ALPMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 91,
    "ALPMIN")), &thisId);
  thisId.fIdentifier = "ALPMAX";
  y->ALPMAX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 92,
    "ALPMAX")), &thisId);
  thisId.fIdentifier = "CAMMIN";
  y->CAMMIN = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 93,
    "CAMMIN")), &thisId);
  thisId.fIdentifier = "CAMMAX";
  y->CAMMAX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 94,
    "CAMMAX")), &thisId);
  thisId.fIdentifier = "LFZO";
  y->LFZO = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 95,
    "LFZO")), &thisId);
  thisId.fIdentifier = "LCX";
  y->LCX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 96,
    "LCX")), &thisId);
  thisId.fIdentifier = "LMUX";
  y->LMUX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 97,
    "LMUX")), &thisId);
  thisId.fIdentifier = "LEX";
  y->LEX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 98,
    "LEX")), &thisId);
  thisId.fIdentifier = "LKX";
  y->LKX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 99,
    "LKX")), &thisId);
  thisId.fIdentifier = "LHX";
  y->LHX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 100,
    "LHX")), &thisId);
  thisId.fIdentifier = "LVX";
  y->LVX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 101,
    "LVX")), &thisId);
  thisId.fIdentifier = "LCY";
  y->LCY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 102,
    "LCY")), &thisId);
  thisId.fIdentifier = "LMUY";
  y->LMUY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 103,
    "LMUY")), &thisId);
  thisId.fIdentifier = "LEY";
  y->LEY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 104,
    "LEY")), &thisId);
  thisId.fIdentifier = "LKY";
  y->LKY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 105,
    "LKY")), &thisId);
  thisId.fIdentifier = "LHY";
  y->LHY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 106,
    "LHY")), &thisId);
  thisId.fIdentifier = "LVY";
  y->LVY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 107,
    "LVY")), &thisId);
  thisId.fIdentifier = "LTR";
  y->LTR = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 108,
    "LTR")), &thisId);
  thisId.fIdentifier = "LRES";
  y->LRES = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 109,
    "LRES")), &thisId);
  thisId.fIdentifier = "LXAL";
  y->LXAL = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 110,
    "LXAL")), &thisId);
  thisId.fIdentifier = "LYKA";
  y->LYKA = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 111,
    "LYKA")), &thisId);
  thisId.fIdentifier = "LVYKA";
  y->LVYKA = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 112,
    "LVYKA")), &thisId);
  thisId.fIdentifier = "LS";
  y->LS = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 113,
    "LS")), &thisId);
  thisId.fIdentifier = "LKYC";
  y->LKYC = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 114,
    "LKYC")), &thisId);
  thisId.fIdentifier = "LKZC";
  y->LKZC = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 115,
    "LKZC")), &thisId);
  thisId.fIdentifier = "LVMX";
  y->LVMX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 116,
    "LVMX")), &thisId);
  thisId.fIdentifier = "LMX";
  y->LMX = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 117,
    "LMX")), &thisId);
  thisId.fIdentifier = "LMY";
  y->LMY = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 118,
    "LMY")), &thisId);
  thisId.fIdentifier = "LMP";
  y->LMP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 119,
    "LMP")), &thisId);
  thisId.fIdentifier = "PCX1";
  y->PCX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 120,
    "PCX1")), &thisId);
  thisId.fIdentifier = "PDX1";
  y->PDX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 121,
    "PDX1")), &thisId);
  thisId.fIdentifier = "PDX2";
  y->PDX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 122,
    "PDX2")), &thisId);
  thisId.fIdentifier = "PDX3";
  y->PDX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 123,
    "PDX3")), &thisId);
  thisId.fIdentifier = "PEX1";
  y->PEX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 124,
    "PEX1")), &thisId);
  thisId.fIdentifier = "PEX2";
  y->PEX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 125,
    "PEX2")), &thisId);
  thisId.fIdentifier = "PEX3";
  y->PEX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 126,
    "PEX3")), &thisId);
  thisId.fIdentifier = "PEX4";
  y->PEX4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 127,
    "PEX4")), &thisId);
  thisId.fIdentifier = "PKX1";
  y->PKX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 128,
    "PKX1")), &thisId);
  thisId.fIdentifier = "PKX2";
  y->PKX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 129,
    "PKX2")), &thisId);
  thisId.fIdentifier = "PKX3";
  y->PKX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 130,
    "PKX3")), &thisId);
  thisId.fIdentifier = "PHX1";
  y->PHX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 131,
    "PHX1")), &thisId);
  thisId.fIdentifier = "PHX2";
  y->PHX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 132,
    "PHX2")), &thisId);
  thisId.fIdentifier = "PVX1";
  y->PVX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 133,
    "PVX1")), &thisId);
  thisId.fIdentifier = "PVX2";
  y->PVX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 134,
    "PVX2")), &thisId);
  thisId.fIdentifier = "PPX1";
  y->PPX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 135,
    "PPX1")), &thisId);
  thisId.fIdentifier = "PPX2";
  y->PPX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 136,
    "PPX2")), &thisId);
  thisId.fIdentifier = "PPX3";
  y->PPX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 137,
    "PPX3")), &thisId);
  thisId.fIdentifier = "PPX4";
  y->PPX4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 138,
    "PPX4")), &thisId);
  thisId.fIdentifier = "RBX1";
  y->RBX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 139,
    "RBX1")), &thisId);
  thisId.fIdentifier = "RBX2";
  y->RBX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 140,
    "RBX2")), &thisId);
  thisId.fIdentifier = "RBX3";
  y->RBX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 141,
    "RBX3")), &thisId);
  thisId.fIdentifier = "RCX1";
  y->RCX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 142,
    "RCX1")), &thisId);
  thisId.fIdentifier = "REX1";
  y->REX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 143,
    "REX1")), &thisId);
  thisId.fIdentifier = "REX2";
  y->REX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 144,
    "REX2")), &thisId);
  thisId.fIdentifier = "RHX1";
  y->RHX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 145,
    "RHX1")), &thisId);
  thisId.fIdentifier = "QSX1";
  y->QSX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 146,
    "QSX1")), &thisId);
  thisId.fIdentifier = "QSX2";
  y->QSX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 147,
    "QSX2")), &thisId);
  thisId.fIdentifier = "QSX3";
  y->QSX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 148,
    "QSX3")), &thisId);
  thisId.fIdentifier = "QSX4";
  y->QSX4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 149,
    "QSX4")), &thisId);
  thisId.fIdentifier = "QSX5";
  y->QSX5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 150,
    "QSX5")), &thisId);
  thisId.fIdentifier = "QSX6";
  y->QSX6 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 151,
    "QSX6")), &thisId);
  thisId.fIdentifier = "QSX7";
  y->QSX7 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 152,
    "QSX7")), &thisId);
  thisId.fIdentifier = "QSX8";
  y->QSX8 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 153,
    "QSX8")), &thisId);
  thisId.fIdentifier = "QSX9";
  y->QSX9 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 154,
    "QSX9")), &thisId);
  thisId.fIdentifier = "QSX10";
  y->QSX10 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 155,
    "QSX10")), &thisId);
  thisId.fIdentifier = "QSX11";
  y->QSX11 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 156,
    "QSX11")), &thisId);
  thisId.fIdentifier = "QSX12";
  y->QSX12 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 157,
    "QSX12")), &thisId);
  thisId.fIdentifier = "QSX13";
  y->QSX13 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 158,
    "QSX13")), &thisId);
  thisId.fIdentifier = "QSX14";
  y->QSX14 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 159,
    "QSX14")), &thisId);
  thisId.fIdentifier = "PPMX1";
  y->PPMX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 160,
    "PPMX1")), &thisId);
  thisId.fIdentifier = "PCY1";
  y->PCY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 161,
    "PCY1")), &thisId);
  thisId.fIdentifier = "PDY1";
  y->PDY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 162,
    "PDY1")), &thisId);
  thisId.fIdentifier = "PDY2";
  y->PDY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 163,
    "PDY2")), &thisId);
  thisId.fIdentifier = "PDY3";
  y->PDY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 164,
    "PDY3")), &thisId);
  thisId.fIdentifier = "PEY1";
  y->PEY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 165,
    "PEY1")), &thisId);
  thisId.fIdentifier = "PEY2";
  y->PEY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 166,
    "PEY2")), &thisId);
  thisId.fIdentifier = "PEY3";
  y->PEY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 167,
    "PEY3")), &thisId);
  thisId.fIdentifier = "PEY4";
  y->PEY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 168,
    "PEY4")), &thisId);
  thisId.fIdentifier = "PEY5";
  y->PEY5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 169,
    "PEY5")), &thisId);
  thisId.fIdentifier = "PKY1";
  y->PKY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 170,
    "PKY1")), &thisId);
  thisId.fIdentifier = "PKY2";
  y->PKY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 171,
    "PKY2")), &thisId);
  thisId.fIdentifier = "PKY3";
  y->PKY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 172,
    "PKY3")), &thisId);
  thisId.fIdentifier = "PKY4";
  y->PKY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 173,
    "PKY4")), &thisId);
  thisId.fIdentifier = "PKY5";
  y->PKY5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 174,
    "PKY5")), &thisId);
  thisId.fIdentifier = "PKY6";
  y->PKY6 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 175,
    "PKY6")), &thisId);
  thisId.fIdentifier = "PKY7";
  y->PKY7 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 176,
    "PKY7")), &thisId);
  thisId.fIdentifier = "PHY1";
  y->PHY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 177,
    "PHY1")), &thisId);
  thisId.fIdentifier = "PHY2";
  y->PHY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 178,
    "PHY2")), &thisId);
  thisId.fIdentifier = "PVY1";
  y->PVY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 179,
    "PVY1")), &thisId);
  thisId.fIdentifier = "PVY2";
  y->PVY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 180,
    "PVY2")), &thisId);
  thisId.fIdentifier = "PVY3";
  y->PVY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 181,
    "PVY3")), &thisId);
  thisId.fIdentifier = "PVY4";
  y->PVY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 182,
    "PVY4")), &thisId);
  thisId.fIdentifier = "PPY1";
  y->PPY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 183,
    "PPY1")), &thisId);
  thisId.fIdentifier = "PPY2";
  y->PPY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 184,
    "PPY2")), &thisId);
  thisId.fIdentifier = "PPY3";
  y->PPY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 185,
    "PPY3")), &thisId);
  thisId.fIdentifier = "PPY4";
  y->PPY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 186,
    "PPY4")), &thisId);
  thisId.fIdentifier = "PPY5";
  y->PPY5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 187,
    "PPY5")), &thisId);
  thisId.fIdentifier = "RBY1";
  y->RBY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 188,
    "RBY1")), &thisId);
  thisId.fIdentifier = "RBY2";
  y->RBY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 189,
    "RBY2")), &thisId);
  thisId.fIdentifier = "RBY3";
  y->RBY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 190,
    "RBY3")), &thisId);
  thisId.fIdentifier = "RBY4";
  y->RBY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 191,
    "RBY4")), &thisId);
  thisId.fIdentifier = "RCY1";
  y->RCY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 192,
    "RCY1")), &thisId);
  thisId.fIdentifier = "REY1";
  y->REY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 193,
    "REY1")), &thisId);
  thisId.fIdentifier = "REY2";
  y->REY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 194,
    "REY2")), &thisId);
  thisId.fIdentifier = "RHY1";
  y->RHY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 195,
    "RHY1")), &thisId);
  thisId.fIdentifier = "RHY2";
  y->RHY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 196,
    "RHY2")), &thisId);
  thisId.fIdentifier = "RVY1";
  y->RVY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 197,
    "RVY1")), &thisId);
  thisId.fIdentifier = "RVY2";
  y->RVY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 198,
    "RVY2")), &thisId);
  thisId.fIdentifier = "RVY3";
  y->RVY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 199,
    "RVY3")), &thisId);
  thisId.fIdentifier = "RVY4";
  y->RVY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 200,
    "RVY4")), &thisId);
  thisId.fIdentifier = "RVY5";
  y->RVY5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 201,
    "RVY5")), &thisId);
  thisId.fIdentifier = "RVY6";
  y->RVY6 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 202,
    "RVY6")), &thisId);
  thisId.fIdentifier = "QSY1";
  y->QSY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 203,
    "QSY1")), &thisId);
  thisId.fIdentifier = "QSY2";
  y->QSY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 204,
    "QSY2")), &thisId);
  thisId.fIdentifier = "QSY3";
  y->QSY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 205,
    "QSY3")), &thisId);
  thisId.fIdentifier = "QSY4";
  y->QSY4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 206,
    "QSY4")), &thisId);
  thisId.fIdentifier = "QSY5";
  y->QSY5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 207,
    "QSY5")), &thisId);
  thisId.fIdentifier = "QSY6";
  y->QSY6 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 208,
    "QSY6")), &thisId);
  thisId.fIdentifier = "QSY7";
  y->QSY7 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 209,
    "QSY7")), &thisId);
  thisId.fIdentifier = "QSY8";
  y->QSY8 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 210,
    "QSY8")), &thisId);
  thisId.fIdentifier = "QBZ1";
  y->QBZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 211,
    "QBZ1")), &thisId);
  thisId.fIdentifier = "QBZ2";
  y->QBZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 212,
    "QBZ2")), &thisId);
  thisId.fIdentifier = "QBZ3";
  y->QBZ3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 213,
    "QBZ3")), &thisId);
  thisId.fIdentifier = "QBZ4";
  y->QBZ4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 214,
    "QBZ4")), &thisId);
  thisId.fIdentifier = "QBZ5";
  y->QBZ5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 215,
    "QBZ5")), &thisId);
  thisId.fIdentifier = "QBZ9";
  y->QBZ9 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 216,
    "QBZ9")), &thisId);
  thisId.fIdentifier = "QBZ10";
  y->QBZ10 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 217,
    "QBZ10")), &thisId);
  thisId.fIdentifier = "QCZ1";
  y->QCZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 218,
    "QCZ1")), &thisId);
  thisId.fIdentifier = "QDZ1";
  y->QDZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 219,
    "QDZ1")), &thisId);
  thisId.fIdentifier = "QDZ2";
  y->QDZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 220,
    "QDZ2")), &thisId);
  thisId.fIdentifier = "QDZ3";
  y->QDZ3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 221,
    "QDZ3")), &thisId);
  thisId.fIdentifier = "QDZ4";
  y->QDZ4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 222,
    "QDZ4")), &thisId);
  thisId.fIdentifier = "QDZ6";
  y->QDZ6 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 223,
    "QDZ6")), &thisId);
  thisId.fIdentifier = "QDZ7";
  y->QDZ7 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 224,
    "QDZ7")), &thisId);
  thisId.fIdentifier = "QDZ8";
  y->QDZ8 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 225,
    "QDZ8")), &thisId);
  thisId.fIdentifier = "QDZ9";
  y->QDZ9 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 226,
    "QDZ9")), &thisId);
  thisId.fIdentifier = "QDZ10";
  y->QDZ10 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 227,
    "QDZ10")), &thisId);
  thisId.fIdentifier = "QDZ11";
  y->QDZ11 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 228,
    "QDZ11")), &thisId);
  thisId.fIdentifier = "QEZ1";
  y->QEZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 229,
    "QEZ1")), &thisId);
  thisId.fIdentifier = "QEZ2";
  y->QEZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 230,
    "QEZ2")), &thisId);
  thisId.fIdentifier = "QEZ3";
  y->QEZ3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 231,
    "QEZ3")), &thisId);
  thisId.fIdentifier = "QEZ4";
  y->QEZ4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 232,
    "QEZ4")), &thisId);
  thisId.fIdentifier = "QEZ5";
  y->QEZ5 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 233,
    "QEZ5")), &thisId);
  thisId.fIdentifier = "QHZ1";
  y->QHZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 234,
    "QHZ1")), &thisId);
  thisId.fIdentifier = "QHZ2";
  y->QHZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 235,
    "QHZ2")), &thisId);
  thisId.fIdentifier = "QHZ3";
  y->QHZ3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 236,
    "QHZ3")), &thisId);
  thisId.fIdentifier = "QHZ4";
  y->QHZ4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 237,
    "QHZ4")), &thisId);
  thisId.fIdentifier = "PPZ1";
  y->PPZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 238,
    "PPZ1")), &thisId);
  thisId.fIdentifier = "PPZ2";
  y->PPZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 239,
    "PPZ2")), &thisId);
  thisId.fIdentifier = "SSZ1";
  y->SSZ1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 240,
    "SSZ1")), &thisId);
  thisId.fIdentifier = "SSZ2";
  y->SSZ2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 241,
    "SSZ2")), &thisId);
  thisId.fIdentifier = "SSZ3";
  y->SSZ3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 242,
    "SSZ3")), &thisId);
  thisId.fIdentifier = "SSZ4";
  y->SSZ4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 243,
    "SSZ4")), &thisId);
  thisId.fIdentifier = "PDXP1";
  y->PDXP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 244,
    "PDXP1")), &thisId);
  thisId.fIdentifier = "PDXP2";
  y->PDXP2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 245,
    "PDXP2")), &thisId);
  thisId.fIdentifier = "PDXP3";
  y->PDXP3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 246,
    "PDXP3")), &thisId);
  thisId.fIdentifier = "PKYP1";
  y->PKYP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 247,
    "PKYP1")), &thisId);
  thisId.fIdentifier = "PDYP1";
  y->PDYP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 248,
    "PDYP1")), &thisId);
  thisId.fIdentifier = "PDYP2";
  y->PDYP2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 249,
    "PDYP2")), &thisId);
  thisId.fIdentifier = "PDYP3";
  y->PDYP3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 250,
    "PDYP3")), &thisId);
  thisId.fIdentifier = "PDYP4";
  y->PDYP4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 251,
    "PDYP4")), &thisId);
  thisId.fIdentifier = "PHYP1";
  y->PHYP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 252,
    "PHYP1")), &thisId);
  thisId.fIdentifier = "PHYP2";
  y->PHYP2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 253,
    "PHYP2")), &thisId);
  thisId.fIdentifier = "PHYP3";
  y->PHYP3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 254,
    "PHYP3")), &thisId);
  thisId.fIdentifier = "PHYP4";
  y->PHYP4 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 255,
    "PHYP4")), &thisId);
  thisId.fIdentifier = "PECP1";
  y->PECP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 256,
    "PECP1")), &thisId);
  thisId.fIdentifier = "PECP2";
  y->PECP2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 257,
    "PECP2")), &thisId);
  thisId.fIdentifier = "QDTP1";
  y->QDTP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 258,
    "QDTP1")), &thisId);
  thisId.fIdentifier = "QCRP1";
  y->QCRP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 259,
    "QCRP1")), &thisId);
  thisId.fIdentifier = "QCRP2";
  y->QCRP2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 260,
    "QCRP2")), &thisId);
  thisId.fIdentifier = "QBRP1";
  y->QBRP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 261,
    "QBRP1")), &thisId);
  thisId.fIdentifier = "QDRP1";
  y->QDRP1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 262,
    "QDRP1")), &thisId);
  thisId.fIdentifier = "FUNCTION_NAME";
  j_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 263,
    "FUNCTION_NAME")), &thisId, y->FUNCTION_NAME);
  thisId.fIdentifier = "SWITCH_INTEG";
  y->SWITCH_INTEG = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 264, "SWITCH_INTEG")), &thisId);
  thisId.fIdentifier = "Q_FCY2";
  y->Q_FCY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    265, "Q_FCY2")), &thisId);
  thisId.fIdentifier = "Q_CAM1";
  y->Q_CAM1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    266, "Q_CAM1")), &thisId);
  thisId.fIdentifier = "Q_CAM2";
  y->Q_CAM2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    267, "Q_CAM2")), &thisId);
  thisId.fIdentifier = "Q_CAM3";
  y->Q_CAM3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    268, "Q_CAM3")), &thisId);
  thisId.fIdentifier = "Q_FYS1";
  y->Q_FYS1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    269, "Q_FYS1")), &thisId);
  thisId.fIdentifier = "Q_FYS2";
  y->Q_FYS2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    270, "Q_FYS2")), &thisId);
  thisId.fIdentifier = "Q_FYS3";
  y->Q_FYS3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    271, "Q_FYS3")), &thisId);
  thisId.fIdentifier = "ENV_C1";
  y->ENV_C1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    272, "ENV_C1")), &thisId);
  thisId.fIdentifier = "ENV_C2";
  y->ENV_C2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    273, "ENV_C2")), &thisId);
  thisId.fIdentifier = "Q_A1";
  y->Q_A1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 274,
    "Q_A1")), &thisId);
  thisId.fIdentifier = "Q_A2";
  y->Q_A2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 275,
    "Q_A2")), &thisId);
  thisId.fIdentifier = "PHY3";
  y->PHY3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 276,
    "PHY3")), &thisId);
  thisId.fIdentifier = "PTX1";
  y->PTX1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 277,
    "PTX1")), &thisId);
  thisId.fIdentifier = "PTX2";
  y->PTX2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 278,
    "PTX2")), &thisId);
  thisId.fIdentifier = "PTX3";
  y->PTX3 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 279,
    "PTX3")), &thisId);
  thisId.fIdentifier = "PTY1";
  y->PTY1 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 280,
    "PTY1")), &thisId);
  thisId.fIdentifier = "PTY2";
  y->PTY2 = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 281,
    "PTY2")), &thisId);
  thisId.fIdentifier = "LSGKP";
  y->LSGKP = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 282,
    "LSGKP")), &thisId);
  thisId.fIdentifier = "LSGAL";
  y->LSGAL = d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 283,
    "LSGAL")), &thisId);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[5]
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[5])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[7]
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[7])
{
  n_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *parameterSource
 *                const char_T *identifier
 *                struct0_T *y
 * Return Type  : void
 */
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *parameterSource, const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(parameterSource), &thisId, y);
  emlrtDestroyArray(&parameterSource);
}

/*
 * Arguments    : const real_T u[30]
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u[30])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 1, 30 };

  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, *(int32_T (*)[2])&iv1[0], 2);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[8]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[8])
{
  o_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[9]
 * Return Type  : void
 */
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[9])
{
  p_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[4]
 * Return Type  : void
 */
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[4])
{
  q_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[6]
 * Return Type  : void
 */
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[6])
{
  r_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                char_T y[39]
 * Return Type  : void
 */
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[39])
{
  s_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *Fz
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Fz, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(Fz), &thisId);
  emlrtDestroyArray(&Fz);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[5]
 * Return Type  : void
 */
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[5])
{
  static const int32_T dims[2] = { 1, 5 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 5);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[7]
 * Return Type  : void
 */
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[7])
{
  static const int32_T dims[2] = { 1, 7 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 7);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[8]
 * Return Type  : void
 */
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[8])
{
  static const int32_T dims[2] = { 1, 8 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 8);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[9]
 * Return Type  : void
 */
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[9])
{
  static const int32_T dims[2] = { 1, 9 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 9);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[4]
 * Return Type  : void
 */
static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[4])
{
  static const int32_T dims[2] = { 1, 4 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 4);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[6]
 * Return Type  : void
 */
static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[6])
{
  static const int32_T dims[2] = { 1, 6 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 6);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                char_T ret[39]
 * Return Type  : void
 */
static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[39])
{
  static const int32_T dims[2] = { 1, 39 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 39);
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray * const prhs[8]
 *                int32_T nlhs
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void mfeval_wrapper_api(const mxArray * const prhs[8], int32_T nlhs, const
  mxArray *plhs[1])
{
  real_T (*outMF)[30];
  struct0_T parameterSource;
  real_T Fz;
  real_T kappa;
  real_T alpha;
  real_T b_gamma;
  real_T phit;
  real_T Vx;
  real_T useMode;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;
  outMF = (real_T (*)[30])mxMalloc(sizeof(real_T [30]));

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "parameterSource",
                   &parameterSource);
  Fz = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Fz");
  kappa = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "kappa");
  alpha = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "alpha");
  b_gamma = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "gamma");
  phit = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "phit");
  Vx = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "Vx");
  useMode = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "useMode");

  /* Invoke the target function */
  mfeval_wrapper(&parameterSource, Fz, kappa, alpha, b_gamma, phit, Vx, useMode,
                 *outMF);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*outMF);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void mfeval_wrapper_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  mfeval_wrapper_xil_terminate();
  mfeval_wrapper_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void mfeval_wrapper_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void mfeval_wrapper_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_mfeval_wrapper_api.c
 *
 * [EOF]
 */
