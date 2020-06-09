/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_mfeval_wrapper_mex.c
 *
 * MATLAB Coder version            : 5.0
 * C/C++ source code generated on  : 09-Jun-2020 11:23:54
 */

/* Include Files */
#include "_coder_mfeval_wrapper_mex.h"
#include "_coder_mfeval_wrapper_api.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void mfeval_wrapper_mexFunction(int32_T nlhs, mxArray *plhs
  [1], int32_T nrhs, const mxArray *prhs[8]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[8]
 * Return Type  : void
 */
void mfeval_wrapper_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[8])
{
  const mxArray *outputs[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        14, "mfeval_wrapper");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 14,
                        "mfeval_wrapper");
  }

  /* Call the function. */
  mfeval_wrapper_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, plhs, outputs);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&mfeval_wrapper_atexit);

  /* Module initialization. */
  mfeval_wrapper_initialize();

  /* Dispatch the entry-point. */
  mfeval_wrapper_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  mfeval_wrapper_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_mfeval_wrapper_mex.c
 *
 * [EOF]
 */
