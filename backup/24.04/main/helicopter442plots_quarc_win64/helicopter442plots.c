/*
 * helicopter442plots.c
 *
 * Code generation for model "helicopter442plots".
 *
 * Model version              : 1.192
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Sun Apr 23 15:01:44 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter442plots.h"
#include "helicopter442plots_private.h"
#include "helicopter442plots_dt.h"

/* Block signals (auto storage) */
B_helicopter442plots_T helicopter442plots_B;

/* Continuous states */
X_helicopter442plots_T helicopter442plots_X;

/* Block states (auto storage) */
DW_helicopter442plots_T helicopter442plots_DW;

/* Real-time model */
RT_MODEL_helicopter442plots_T helicopter442plots_M_;
RT_MODEL_helicopter442plots_T *const helicopter442plots_M =
  &helicopter442plots_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter442plots_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter442plots_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWorkspace[2];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[12];
  real_T *lastU;
  real_T rtb_Derivative;
  int32_T i;
  real_T tmp[6];
  int32_T i_0;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* set solver stop time */
    if (!(helicopter442plots_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter442plots_M->solverInfo,
                            ((helicopter442plots_M->Timing.clockTickH0 + 1) *
        helicopter442plots_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter442plots_M->solverInfo,
                            ((helicopter442plots_M->Timing.clockTick0 + 1) *
        helicopter442plots_M->Timing.stepSize0 +
        helicopter442plots_M->Timing.clockTickH0 *
        helicopter442plots_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter442plots_M)) {
    helicopter442plots_M->Timing.t[0] = rtsiGetT
      (&helicopter442plots_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter442plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter442plots_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter442plots_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter442plots_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter442plots_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter442plots_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter442plots_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter442plots_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter442plots_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter442plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[120]) {
      currTimeIndex = 119;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter442plots_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&helicopter442plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 121;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&helicopter442plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 121;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter442plots_B.FromWorkspace1[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 121;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter442plots_B.TravelCounttorad =
      helicopter442plots_P.TravelCounttorad_Gain * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helicopter442plots_B.Gain = helicopter442plots_P.Gain_Gain *
      helicopter442plots_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    helicopter442plots_B.Sum1 = helicopter442plots_P.travel_offsetdeg_Value +
      helicopter442plots_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter442plots_P.TravelTransferFcn_C *
    helicopter442plots_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter442plots_P.TravelTransferFcn_D *
    helicopter442plots_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helicopter442plots_B.Gain_d = helicopter442plots_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter442plots_B.PitchCounttorad =
      helicopter442plots_P.PitchCounttorad_Gain * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helicopter442plots_B.Gain_i = helicopter442plots_P.Gain_Gain_a *
      helicopter442plots_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter442plots_P.PitchTransferFcn_C *
    helicopter442plots_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter442plots_P.PitchTransferFcn_D *
    helicopter442plots_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  helicopter442plots_B.Gain_b = helicopter442plots_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter442plots_B.ElevationCounttorad =
      helicopter442plots_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helicopter442plots_B.Gain_e = helicopter442plots_P.Gain_Gain_lv *
      helicopter442plots_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter442plots_B.Sum = helicopter442plots_B.Gain_e +
      helicopter442plots_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter442plots_P.ElevationTransferFcn_C *
    helicopter442plots_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter442plots_P.ElevationTransferFcn_D *
    helicopter442plots_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helicopter442plots_B.Gain_dg = helicopter442plots_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter442plots_B.Gain1[0] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Sum1;
  helicopter442plots_B.Gain1[1] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Gain_d;
  helicopter442plots_B.Gain1[2] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Gain_i;
  helicopter442plots_B.Gain1[3] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Gain_b;
  helicopter442plots_B.Gain1[4] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Sum;
  helicopter442plots_B.Gain1[5] = helicopter442plots_P.Gain1_Gain *
    helicopter442plots_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i] = helicopter442plots_B.FromWorkspace1[i];
    }

    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i + 6] = helicopter442plots_B.Gain1[i];
    }

    /* End of SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter442plots_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter442plots_DW.ToFile_IWORK.Count*13)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter442plots_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[13];
          helicopter442plots_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter442plots_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          u[7] = rtb_TmpSignalConversionAtToFile[6];
          u[8] = rtb_TmpSignalConversionAtToFile[7];
          u[9] = rtb_TmpSignalConversionAtToFile[8];
          u[10] = rtb_TmpSignalConversionAtToFile[9];
          u[11] = rtb_TmpSignalConversionAtToFile[10];
          u[12] = rtb_TmpSignalConversionAtToFile[11];
          if (fwrite(u, sizeof(real_T), 13, fp) != 13) {
            rtmSetErrorStatus(helicopter442plots_M,
                              "Error writing to MAT-file 442.mat");
            return;
          }

          if (((++helicopter442plots_DW.ToFile_IWORK.Count)*13)+1 >= 100000000)
          {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file 442.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter442plots_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter442plots_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter442plots_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter442plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[120]) {
      currTimeIndex = 119;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter442plots_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 121;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 121;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWorkspace[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 121;
          }
        }
      }
    }
  }

  /* Sum: '<S5>/Sum1' incorporates:
   *  Gain: '<S5>/Gain'
   */
  for (i = 0; i < 6; i++) {
    tmp[i] = helicopter442plots_B.Gain1[i] -
      helicopter442plots_B.FromWorkspace1[i];
  }

  /* End of Sum: '<S5>/Sum1' */

  /* Sum: '<S5>/Sum' incorporates:
   *  Gain: '<S5>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Derivative = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Derivative += helicopter442plots_P.K[(i_0 << 1) + i] * tmp[i_0];
    }

    helicopter442plots_B.Sum_g[i] = rtb_FromWorkspace[i] - rtb_Derivative;
  }

  /* End of Sum: '<S5>/Sum' */

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<S6>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helicopter442plots_B.Sum_k = ((helicopter442plots_B.Sum_g[0] -
    helicopter442plots_B.Gain1[2]) * helicopter442plots_P.K_pp -
    helicopter442plots_P.K_pd * helicopter442plots_B.Gain1[3]) +
    helicopter442plots_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter442plots_X.Integrator_CSTATE >=
      helicopter442plots_P.Integrator_UpperSat ) {
    helicopter442plots_X.Integrator_CSTATE =
      helicopter442plots_P.Integrator_UpperSat;
  } else if (helicopter442plots_X.Integrator_CSTATE <=
             (helicopter442plots_P.Integrator_LowerSat) ) {
    helicopter442plots_X.Integrator_CSTATE =
      (helicopter442plots_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter442plots_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' */
  rtb_Derivative = helicopter442plots_B.Sum_g[1] - helicopter442plots_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter442plots_B.Sum2 = ((helicopter442plots_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter442plots_P.K_ed * helicopter442plots_B.Gain1[5]) +
    helicopter442plots_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace2' */
  {
    real_T *pDataValues = (real_T *)
      helicopter442plots_DW.FromWorkspace2_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter442plots_DW.FromWorkspace2_PWORK.TimePtr;
    int_T currTimeIndex = helicopter442plots_DW.FromWorkspace2_IWORK.PrevIndex;
    real_T t = helicopter442plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[120]) {
      currTimeIndex = 119;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter442plots_DW.FromWorkspace2_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helicopter442plots_B.FromWorkspace2 = pDataValues[currTimeIndex];
        } else {
          helicopter442plots_B.FromWorkspace2 = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helicopter442plots_B.FromWorkspace2 = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 121;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter442plots_B.Sum2 - helicopter442plots_B.Sum_k) *
    helicopter442plots_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter442plots_B.K_ei = helicopter442plots_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter442plots_DW.TimeStampA >= helicopter442plots_M->Timing.t[0]) &&
      (helicopter442plots_DW.TimeStampB >= helicopter442plots_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter442plots_DW.TimeStampA;
    lastU = &helicopter442plots_DW.LastUAtTimeA;
    if (helicopter442plots_DW.TimeStampA < helicopter442plots_DW.TimeStampB) {
      if (helicopter442plots_DW.TimeStampB < helicopter442plots_M->Timing.t[0])
      {
        rtb_Derivative = helicopter442plots_DW.TimeStampB;
        lastU = &helicopter442plots_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter442plots_DW.TimeStampA >= helicopter442plots_M->Timing.t[0])
      {
        rtb_Derivative = helicopter442plots_DW.TimeStampB;
        lastU = &helicopter442plots_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter442plots_B.PitchCounttorad - *lastU) /
      (helicopter442plots_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helicopter442plots_B.Gain_l = helicopter442plots_P.Gain_Gain_a1 *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter442plots_P.BackmotorSaturation_UpperSat) {
    helicopter442plots_B.BackmotorSaturation =
      helicopter442plots_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter442plots_P.BackmotorSaturation_LowerSat) {
    helicopter442plots_B.BackmotorSaturation =
      helicopter442plots_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter442plots_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helicopter442plots_B.Sum_k + helicopter442plots_B.Sum2) *
    helicopter442plots_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter442plots_P.FrontmotorSaturation_UpperSat) {
    helicopter442plots_B.FrontmotorSaturation =
      helicopter442plots_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter442plots_P.FrontmotorSaturation_LowerSat)
  {
    helicopter442plots_B.FrontmotorSaturation =
      helicopter442plots_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter442plots_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter442plots/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter442plots_DW.HILWriteAnalog_Buffer[0] =
        helicopter442plots_B.FrontmotorSaturation;
      helicopter442plots_DW.HILWriteAnalog_Buffer[1] =
        helicopter442plots_B.BackmotorSaturation;
      result = hil_write_analog(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILWriteAnalog_channels, 2,
        &helicopter442plots_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter442plots_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter442plots_DW.TimeStampA == (rtInf)) {
    helicopter442plots_DW.TimeStampA = helicopter442plots_M->Timing.t[0];
    lastU = &helicopter442plots_DW.LastUAtTimeA;
  } else if (helicopter442plots_DW.TimeStampB == (rtInf)) {
    helicopter442plots_DW.TimeStampB = helicopter442plots_M->Timing.t[0];
    lastU = &helicopter442plots_DW.LastUAtTimeB;
  } else if (helicopter442plots_DW.TimeStampA < helicopter442plots_DW.TimeStampB)
  {
    helicopter442plots_DW.TimeStampA = helicopter442plots_M->Timing.t[0];
    lastU = &helicopter442plots_DW.LastUAtTimeA;
  } else {
    helicopter442plots_DW.TimeStampB = helicopter442plots_M->Timing.t[0];
    lastU = &helicopter442plots_DW.LastUAtTimeB;
  }

  *lastU = helicopter442plots_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter442plots_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter442plots_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter442plots_M->Timing.clockTick0)) {
    ++helicopter442plots_M->Timing.clockTickH0;
  }

  helicopter442plots_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter442plots_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter442plots_M->Timing.clockTick1)) {
      ++helicopter442plots_M->Timing.clockTickH1;
    }

    helicopter442plots_M->Timing.t[1] = helicopter442plots_M->Timing.clockTick1 *
      helicopter442plots_M->Timing.stepSize1 +
      helicopter442plots_M->Timing.clockTickH1 *
      helicopter442plots_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter442plots_derivatives(void)
{
  XDot_helicopter442plots_T *_rtXdot;
  _rtXdot = ((XDot_helicopter442plots_T *)
             helicopter442plots_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter442plots_P.TravelTransferFcn_A *
    helicopter442plots_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter442plots_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter442plots_P.PitchTransferFcn_A *
    helicopter442plots_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter442plots_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter442plots_P.ElevationTransferFcn_A *
    helicopter442plots_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter442plots_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter442plots_X.Integrator_CSTATE <=
            (helicopter442plots_P.Integrator_LowerSat) );
    usat = ( helicopter442plots_X.Integrator_CSTATE >=
            helicopter442plots_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter442plots_B.K_ei > 0)) ||
        (usat && (helicopter442plots_B.K_ei < 0)) ) {
      ((XDot_helicopter442plots_T *) helicopter442plots_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter442plots_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter442plots_T *) helicopter442plots_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter442plots_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter442plots/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter442plots_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter442plots_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter442plots_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
      return;
    }

    if ((helicopter442plots_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_analog_inpu_m && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter442plots_DW.HILInitialize_AIMinimums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter442plots_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter442plots_DW.HILInitialize_AIMaximums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter442plots_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_analog_input_chan, 8U,
         &helicopter442plots_DW.HILInitialize_AIMinimums[0],
         &helicopter442plots_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_analog_outp_b && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter442plots_DW.HILInitialize_AOMinimums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter442plots_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter442plots_DW.HILInitialize_AOMaximums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter442plots_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter442plots_DW.HILInitialize_AOMinimums[0],
         &helicopter442plots_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_analog_outp_j && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter442plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter442plots_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILInitialize_analog_output_cha, 8U,
        &helicopter442plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter442plots_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter442plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter442plots_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter442plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_encoder_par_m && is_switching))
    {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter442plots_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter442plots_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter442plots_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_encoder_cou_k && is_switching))
    {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter442plots_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter442plots_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILInitialize_encoder_channels, 8U,
        &helicopter442plots_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_pwm_params__f && is_switching))
    {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter442plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter442plots_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter442plots_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter442plots_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter442plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter442plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter442plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter442plots_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter442plots_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter442plots_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter442plots_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter442plots_DW.HILInitialize_Card,
          &helicopter442plots_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter442plots_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter442plots_DW.HILInitialize_Card,
          &helicopter442plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter442plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter442plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter442plots_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter442plots_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            helicopter442plots_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter442plots_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            helicopter442plots_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_pwm_channels, 8U,
         (t_pwm_configuration *)
         &helicopter442plots_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *) &helicopter442plots_DW.HILInitialize_POAlignValues
         [0],
         (t_pwm_polarity *) &helicopter442plots_DW.HILInitialize_POPolarityVals
         [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter442plots_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter442plots_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter442plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter442plots_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter442plots_DW.HILInitialize_POSortedFreqs[0],
        &helicopter442plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter442plots_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_pwm_outputs_g && is_switching))
    {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter442plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter442plots_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter442plots_DW.HILInitialize_Card,
        helicopter442plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter442plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter442plots_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter442plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter442plots_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter442plots_DW.HILInitialize_Card,
         helicopter442plots_P.HILInitialize_pwm_channels, 8U,
         &helicopter442plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter442plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter442plots_DW.HILInitialize_Card,
       helicopter442plots_P.HILReadEncoderTimebase_samples_,
       helicopter442plots_P.HILReadEncoderTimebase_channels, 3,
       &helicopter442plots_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      2.6261999647428622, 2.0911635140141875, 1.6622992678721005,
      1.5877238717310451, 1.573471716467016, 1.5511246149350129,
      1.515576787795625, 1.4653752785251928, 1.4000871274543616,
      1.3205254992527797, 1.22869887857918, 1.1274866599125919,
      1.0202120144514262, 0.91027700239572329, 0.80086052722204748,
      0.69472172069789, 0.59411377997758485, 0.50078315632898185,
      0.41601808099384097, 0.34071906513697603, 0.27547299875770154,
      0.22062379132718068, 0.1763342334929503, 0.14051524038265983,
      0.11081957446577165, 0.085297222632665759, 0.062513209815632673,
      0.041478942406035807, 0.021541047392226564, 0.0022798640280272558,
      -0.016568589171363848, -0.035167744463503681, -0.053617522932316308,
      -0.071978399512719241, -0.090286670504554342, -0.10856396817062969,
      -0.12682310836868357, -0.14507164478288717, -0.163314009669106,
      -0.18155279351904169, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1078830433320031,
      -0.035551938161722259, -0.016874558104801, -0.045491798709801856,
      -0.089240050847105887, -0.14224357596649606, -0.20083328084420315,
      -0.2611606303727394, -0.31824785902687192, -0.36730661125129171,
      -0.40484888197332264, -0.42909858216701835, -0.43974004823182777,
      -0.43766590069477795, -0.4245552260969141, -0.40243176288186716,
      -0.3733224945935022, -0.3390603013408206, -0.30119606342795835,
      -0.26098426551740134, -0.21939682972189958, -0.17715823133677872,
      -0.14327597244125881, -0.11878266366745857, -0.10208940733255685,
      -0.091136051268231086, -0.084137069638456963, -0.079751580055369933,
      -0.077044733457033346, -0.075393812797515433, -0.074396621168279276,
      -0.073799113875490888, -0.073443506321396984, -0.0732330839679359,
      -0.073109190664188736, -0.073036560791999577, -0.072994145656977386,
      -0.072969459544798207, -0.072955135399750293, -0.072946846389123132, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -0.0022554010138257242, 0.10352680562479485,
      0.21730495019445994, 0.30957294386705281, 0.37461314579700611,
      0.41408910437420127, 0.42637025230886483, 0.40347058696813787,
      0.34672849873422618, 0.26533442879002889, 0.17138761931152721,
      0.075209817917711988, -0.014659282624286016, -0.092661240785411755,
      -0.15636018854090888, -0.20573319075669327, -0.24215209638388366,
      -0.26760997233684009, -0.28420162981431213, -0.29392411299556787,
      -0.29852628148562516, -0.23946686545302248, -0.17310935184891441,
      -0.11798156023427767, -0.07741413732341515, -0.049466129084870088,
      -0.030994965453042242, -0.019130957950166318, -0.011668076695620418,
      -0.0070477695852658763, -0.0042229533433716625, -0.0025132983781613947,
      -0.0014871848296109209, -0.00087563055317940134, -0.000513320197540521,
      -0.00029977397495733587, -0.00017447201269744247, -0.00010123758654236517,
      -5.858356135560049E-5, -3.3817547562146827E-5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.41530526505749382, 0.45405456352359957, 0.369034126731209,
      0.26016021579332083, 0.15790382904947439, 0.049124591987844239,
      -0.091598660940133331, -0.22696835275631846, -0.32557627974263209,
      -0.37578723791078772, -0.38471120557512978, -0.35947640216803134,
      -0.31200783264452081, -0.25479579102200334, -0.19749200886313095,
      -0.14567562250874866, -0.1018315038118231, -0.066366629909904717,
      -0.038889932725000249, -0.018408673960258359, 0.23623766413040947,
      0.2654300544163668, 0.22051116645863128, 0.1622696916435051,
      0.11179203295412675, 0.07388465452739984, 0.047456030011488314,
      0.029851525018109146, 0.018481228441443783, 0.01129926496759661,
      0.0068386198608344955, 0.0041044541942076443, 0.002446217105731654,
      0.0014492414225604805, 0.00085418489033062684, 0.00050120784904000978,
      0.00029293770461942862, 0.00017061610074675258, 9.9064055173959936E-5,
      5.7359090215455356E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15351669481296898,
      0.069192994718704048, 0.067998765196741723, 0.076615701361009925,
      0.085300282373037759, 0.093713949998177545, 0.10174587200612085,
      0.10928217747852303, 0.11636513990158449, 0.12289924113462876,
      0.12889960533412081, 0.13437907925820816, 0.13937772871387039,
      0.14391851781546014, 0.14802974488096349, 0.15174110004662125,
      0.15508067818733867, 0.1580725749367885, 0.16074440582127647,
      0.16312768993449445, 0.1652057513704559, 0.16500987140308052,
      0.16307610250069321, 0.15984197365045114, 0.15566278920811097,
      0.15082534549668458, 0.14555950621058666, 0.14004795446631718,
      0.13443439637068055, 0.12883045086654882, 0.12332142611714036,
      0.11797115314751633, 0.11282602217751406, 0.10791834545035368,
      0.10326915186806793, 0.098890503129342069, 0.094787424423530678,
      0.0909602634071056, 0.087419122863070936, 0.084170721660241546, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.032924301503279206, 0.034241242822684856,
      0.034988593108849657, 0.034447288052041752, 0.033371595922343983,
      0.031831982285359439, 0.029817945983796463, 0.027980557039973177,
      0.0257757792685031, 0.023638576141137242, 0.021554350774452665,
      0.019630637250321527, 0.017798884902905224, 0.016080428349667802,
      0.014480837575127596, 0.012993731390371585, 0.011603110252999059,
      0.010323058349618901, 0.0091693295676552253, 0.0079489866140057524,
      -0.0011488713929645785, -0.0081017150235909476, -0.013304106230829522,
      -0.017085033634500438, -0.019718576008048515, -0.021432502620065949,
      -0.022415568310262073, -0.022823708203926273, -0.022785293346664976,
      -0.022405585238258321, -0.021770507421704176, -0.020949835261532861,
      -0.01999989043309431, -0.018965814104768112, -0.017883482038920283,
      -0.016781115773815563, -0.015680663195606293, -0.014601185262726826,
      -0.01363709976541168, -0.013521122312665916, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helicopter442plots_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter442plots_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter442plots_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "442.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter442plots_M, "Error creating .mat file 442.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,13,0,"x_state")) {
      rtmSetErrorStatus(helicopter442plots_M,
                        "Error writing mat file header to file 442.mat");
      return;
    }

    helicopter442plots_DW.ToFile_IWORK.Count = 0;
    helicopter442plots_DW.ToFile_IWORK.Decimation = -1;
    helicopter442plots_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.51270335137808509, 0.50912591265272644, 0.50336431114883806,
      0.49294183913233969, 0.47239757996907306, 0.41576649175309016,
      0.294939573120584, 0.15747088590074002, 0.029545346831469588,
      -0.077011760456090136, -0.163224190697314, -0.22491519749162589,
      -0.26560523446218715, -0.29070250207805282, -0.30502226518162356,
      -0.311824906593301, -0.31346645194305006, -0.31151441678971387,
      -0.30742882386030806, -0.30212716152759422, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.5235987755982987, 0.163306183711081,
      0.15715084051896197, 0.14962694294692852, 0.14784310454974225,
      0.14482096270240313, 0.13950512213427455, 0.14232478850304547,
      0.13664344830899067, 0.13599247915689444, 0.13482686640932565,
      0.13505934072680176, 0.13431598506355127, 0.13380089995105804,
      0.13337018091921438, 0.13288321300545705, 0.13219029409054472,
      0.13173688400687286, 0.13162995033687575, 0.12861406742231732, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter442plots_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter442plots_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter442plots_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace2' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.51270335137808509, 0.50912591265272644, 0.50336431114883806,
      0.49294183913233969, 0.47239757996907306, 0.41576649175309016,
      0.294939573120584, 0.15747088590074002, 0.029545346831469588,
      -0.077011760456090136, -0.163224190697314, -0.22491519749162589,
      -0.26560523446218715, -0.29070250207805282, -0.30502226518162356,
      -0.311824906593301, -0.31346645194305006, -0.31151441678971387,
      -0.30742882386030806, -0.30212716152759422, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter442plots_DW.FromWorkspace2_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter442plots_DW.FromWorkspace2_PWORK.DataPtr = (void *) pDataValues0;
    helicopter442plots_DW.FromWorkspace2_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter442plots_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter442plots_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter442plots_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter442plots_X.Integrator_CSTATE = helicopter442plots_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter442plots_DW.TimeStampA = (rtInf);
  helicopter442plots_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter442plots_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter442plots/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter442plots_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter442plots_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter442plots_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_analog_outp_c && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter442plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter442plots_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter442plots_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter442plots_P.HILInitialize_set_pwm_outputs_p && is_switching))
    {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter442plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter442plots_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter442plots_DW.HILInitialize_Card
                         , helicopter442plots_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter442plots_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter442plots_DW.HILInitialize_AOVoltages[0]
                         , &helicopter442plots_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter442plots_DW.HILInitialize_Card,
             helicopter442plots_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter442plots_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter442plots_DW.HILInitialize_Card,
            helicopter442plots_P.HILInitialize_pwm_channels,
            num_final_pwm_outputs,
            &helicopter442plots_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter442plots_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter442plots_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter442plots_DW.HILInitialize_Card);
    hil_close(helicopter442plots_DW.HILInitialize_Card);
    helicopter442plots_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter442plots_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "442.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter442plots_M, "Error closing MAT-file 442.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter442plots_M,
                          "Error reopening MAT-file 442.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 13,
           helicopter442plots_DW.ToFile_IWORK.Count, "x_state")) {
        rtmSetErrorStatus(helicopter442plots_M,
                          "Error writing header for x_state to MAT-file 442.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter442plots_M, "Error closing MAT-file 442.mat");
        return;
      }

      helicopter442plots_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter442plots_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter442plots_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter442plots_initialize();
}

void MdlTerminate(void)
{
  helicopter442plots_terminate();
}

/* Registration function */
RT_MODEL_helicopter442plots_T *helicopter442plots(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter442plots_P.Integrator_UpperSat = rtInf;
  helicopter442plots_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter442plots_M, 0,
                sizeof(RT_MODEL_helicopter442plots_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter442plots_M->solverInfo,
                          &helicopter442plots_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter442plots_M->solverInfo, &rtmGetTPtr
                (helicopter442plots_M));
    rtsiSetStepSizePtr(&helicopter442plots_M->solverInfo,
                       &helicopter442plots_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter442plots_M->solverInfo,
                 &helicopter442plots_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter442plots_M->solverInfo, (real_T **)
                         &helicopter442plots_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter442plots_M->solverInfo,
      &helicopter442plots_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter442plots_M->solverInfo, (&rtmGetErrorStatus
                           (helicopter442plots_M)));
    rtsiSetRTModelPtr(&helicopter442plots_M->solverInfo, helicopter442plots_M);
  }

  rtsiSetSimTimeStep(&helicopter442plots_M->solverInfo, MAJOR_TIME_STEP);
  helicopter442plots_M->ModelData.intgData.f[0] =
    helicopter442plots_M->ModelData.odeF[0];
  helicopter442plots_M->ModelData.contStates = ((real_T *) &helicopter442plots_X);
  rtsiSetSolverData(&helicopter442plots_M->solverInfo, (void *)
                    &helicopter442plots_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter442plots_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter442plots_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter442plots_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter442plots_M->Timing.sampleTimes =
      (&helicopter442plots_M->Timing.sampleTimesArray[0]);
    helicopter442plots_M->Timing.offsetTimes =
      (&helicopter442plots_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter442plots_M->Timing.sampleTimes[0] = (0.0);
    helicopter442plots_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter442plots_M->Timing.offsetTimes[0] = (0.0);
    helicopter442plots_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter442plots_M, &helicopter442plots_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter442plots_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter442plots_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter442plots_M, -1);
  helicopter442plots_M->Timing.stepSize0 = 0.002;
  helicopter442plots_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter442plots_M->Sizes.checksums[0] = (1559780280U);
  helicopter442plots_M->Sizes.checksums[1] = (2885554622U);
  helicopter442plots_M->Sizes.checksums[2] = (1423887505U);
  helicopter442plots_M->Sizes.checksums[3] = (2261873220U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter442plots_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter442plots_M->extModeInfo,
      &helicopter442plots_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter442plots_M->extModeInfo,
                        helicopter442plots_M->Sizes.checksums);
    rteiSetTPtr(helicopter442plots_M->extModeInfo, rtmGetTPtr
                (helicopter442plots_M));
  }

  helicopter442plots_M->solverInfoPtr = (&helicopter442plots_M->solverInfo);
  helicopter442plots_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter442plots_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter442plots_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter442plots_M->ModelData.blockIO = ((void *) &helicopter442plots_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter442plots_B.FromWorkspace1[i] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      helicopter442plots_B.Gain1[i] = 0.0;
    }

    helicopter442plots_B.TravelCounttorad = 0.0;
    helicopter442plots_B.Gain = 0.0;
    helicopter442plots_B.Sum1 = 0.0;
    helicopter442plots_B.Gain_d = 0.0;
    helicopter442plots_B.PitchCounttorad = 0.0;
    helicopter442plots_B.Gain_i = 0.0;
    helicopter442plots_B.Gain_b = 0.0;
    helicopter442plots_B.ElevationCounttorad = 0.0;
    helicopter442plots_B.Gain_e = 0.0;
    helicopter442plots_B.Sum = 0.0;
    helicopter442plots_B.Gain_dg = 0.0;
    helicopter442plots_B.Sum_g[0] = 0.0;
    helicopter442plots_B.Sum_g[1] = 0.0;
    helicopter442plots_B.Sum_k = 0.0;
    helicopter442plots_B.Sum2 = 0.0;
    helicopter442plots_B.FromWorkspace2 = 0.0;
    helicopter442plots_B.K_ei = 0.0;
    helicopter442plots_B.Gain_l = 0.0;
    helicopter442plots_B.BackmotorSaturation = 0.0;
    helicopter442plots_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter442plots_M->ModelData.defaultParam = ((real_T *)
    &helicopter442plots_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter442plots_X;
    helicopter442plots_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter442plots_X, 0,
                  sizeof(X_helicopter442plots_T));
  }

  /* states (dwork) */
  helicopter442plots_M->ModelData.dwork = ((void *) &helicopter442plots_DW);
  (void) memset((void *)&helicopter442plots_DW, 0,
                sizeof(DW_helicopter442plots_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter442plots_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter442plots_DW.TimeStampA = 0.0;
  helicopter442plots_DW.LastUAtTimeA = 0.0;
  helicopter442plots_DW.TimeStampB = 0.0;
  helicopter442plots_DW.LastUAtTimeB = 0.0;
  helicopter442plots_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter442plots_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter442plots_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter442plots_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter442plots_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter442plots_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter442plots_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter442plots_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter442plots_M->Sizes.numBlocks = (62);/* Number of blocks */
  helicopter442plots_M->Sizes.numBlockIO = (21);/* Number of block outputs */
  helicopter442plots_M->Sizes.numBlockPrms = (153);/* Sum of parameter "widths" */
  return helicopter442plots_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
