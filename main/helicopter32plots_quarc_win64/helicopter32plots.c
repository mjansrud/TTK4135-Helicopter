/*
 * helicopter32plots.c
 *
 * Code generation for model "helicopter32plots".
 *
 * Model version              : 1.180
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Sun Apr 23 13:12:02 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter32plots.h"
#include "helicopter32plots_private.h"
#include "helicopter32plots_dt.h"

/* Block signals (auto storage) */
B_helicopter32plots_T helicopter32plots_B;

/* Continuous states */
X_helicopter32plots_T helicopter32plots_X;

/* Block states (auto storage) */
DW_helicopter32plots_T helicopter32plots_DW;

/* Real-time model */
RT_MODEL_helicopter32plots_T helicopter32plots_M_;
RT_MODEL_helicopter32plots_T *const helicopter32plots_M = &helicopter32plots_M_;

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
  helicopter32plots_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter32plots_output(void)
{
  /* local block i/o variables */
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[10];
  real_T *lastU;
  real_T rtb_Derivative;
  int32_T i;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* set solver stop time */
    if (!(helicopter32plots_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter32plots_M->solverInfo,
                            ((helicopter32plots_M->Timing.clockTickH0 + 1) *
        helicopter32plots_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter32plots_M->solverInfo,
                            ((helicopter32plots_M->Timing.clockTick0 + 1) *
        helicopter32plots_M->Timing.stepSize0 +
        helicopter32plots_M->Timing.clockTickH0 *
        helicopter32plots_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter32plots_M)) {
    helicopter32plots_M->Timing.t[0] = rtsiGetT(&helicopter32plots_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter32plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter32plots_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter32plots_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter32plots_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter32plots_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter32plots_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter32plots_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter32plots_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter32plots_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter32plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helicopter32plots_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter32plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter32plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter32plots_B.FromWorkspace1[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter32plots_B.TravelCounttorad =
      helicopter32plots_P.TravelCounttorad_Gain * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helicopter32plots_B.Gain = helicopter32plots_P.Gain_Gain *
      helicopter32plots_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    helicopter32plots_B.Sum1 = helicopter32plots_P.travel_offsetdeg_Value +
      helicopter32plots_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter32plots_P.TravelTransferFcn_C *
    helicopter32plots_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter32plots_P.TravelTransferFcn_D *
    helicopter32plots_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helicopter32plots_B.Gain_d = helicopter32plots_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter32plots_B.PitchCounttorad =
      helicopter32plots_P.PitchCounttorad_Gain * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helicopter32plots_B.Gain_i = helicopter32plots_P.Gain_Gain_a *
      helicopter32plots_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter32plots_P.PitchTransferFcn_C *
    helicopter32plots_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter32plots_P.PitchTransferFcn_D *
    helicopter32plots_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  helicopter32plots_B.Gain_b = helicopter32plots_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter32plots_B.ElevationCounttorad =
      helicopter32plots_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helicopter32plots_B.Gain_e = helicopter32plots_P.Gain_Gain_lv *
      helicopter32plots_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter32plots_B.Sum = helicopter32plots_B.Gain_e +
      helicopter32plots_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter32plots_P.ElevationTransferFcn_C *
    helicopter32plots_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter32plots_P.ElevationTransferFcn_D *
    helicopter32plots_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helicopter32plots_B.Gain_dg = helicopter32plots_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter32plots_B.Gain1[0] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Sum1;
  helicopter32plots_B.Gain1[1] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Gain_d;
  helicopter32plots_B.Gain1[2] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Gain_i;
  helicopter32plots_B.Gain1[3] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Gain_b;
  helicopter32plots_B.Gain1[4] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Sum;
  helicopter32plots_B.Gain1[5] = helicopter32plots_P.Gain1_Gain *
    helicopter32plots_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    rtb_TmpSignalConversionAtToFile[0] = helicopter32plots_B.FromWorkspace1[0];
    rtb_TmpSignalConversionAtToFile[1] = helicopter32plots_B.FromWorkspace1[1];
    rtb_TmpSignalConversionAtToFile[2] = helicopter32plots_B.FromWorkspace1[2];
    rtb_TmpSignalConversionAtToFile[3] = helicopter32plots_B.FromWorkspace1[3];
    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i + 4] = helicopter32plots_B.Gain1[i];
    }

    /* End of SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter32plots_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter32plots_DW.ToFile_IWORK.Count*11)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter32plots_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[11];
          helicopter32plots_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter32plots_M->Timing.t[1];
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
          if (fwrite(u, sizeof(real_T), 11, fp) != 11) {
            rtmSetErrorStatus(helicopter32plots_M,
                              "Error writing to MAT-file 32.mat");
            return;
          }

          if (((++helicopter32plots_DW.ToFile_IWORK.Count)*11)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file 32.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter32plots_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter32plots_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter32plots_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter32plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helicopter32plots_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helicopter32plots_B.FromWorkspace = pDataValues[currTimeIndex];
        } else {
          helicopter32plots_B.FromWorkspace = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helicopter32plots_B.FromWorkspace = (real_T) rtInterpolate(d1, d2, f1,
          f2);
        pDataValues += 141;
      }
    }
  }

  /* Gain: '<S5>/Gain' incorporates:
   *  Sum: '<S5>/Sum1'
   */
  rtb_Backgain = (((helicopter32plots_B.Gain1[0] -
                    helicopter32plots_B.FromWorkspace1[0]) *
                   helicopter32plots_P.K[0] + (helicopter32plots_B.Gain1[1] -
    helicopter32plots_B.FromWorkspace1[1]) * helicopter32plots_P.K[1]) +
                  (helicopter32plots_B.Gain1[2] -
                   helicopter32plots_B.FromWorkspace1[2]) *
                  helicopter32plots_P.K[2]) + (helicopter32plots_B.Gain1[3] -
    helicopter32plots_B.FromWorkspace1[3]) * helicopter32plots_P.K[3];

  /* Sum: '<S5>/Sum' */
  helicopter32plots_B.Sum_g = helicopter32plots_B.FromWorkspace - rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo WorkspaceInport1' */
    helicopter32plots_B.TmpSignalConversionAtToWorkspac[0] =
      helicopter32plots_B.FromWorkspace;
    helicopter32plots_B.TmpSignalConversionAtToWorkspac[1] =
      helicopter32plots_B.Sum_g;
  }

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<S6>/Vd_bias'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helicopter32plots_B.Sum_k = ((helicopter32plots_B.Sum_g -
    helicopter32plots_B.Gain1[2]) * helicopter32plots_P.K_pp -
    helicopter32plots_P.K_pd * helicopter32plots_B.Gain1[3]) +
    helicopter32plots_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter32plots_X.Integrator_CSTATE >=
      helicopter32plots_P.Integrator_UpperSat ) {
    helicopter32plots_X.Integrator_CSTATE =
      helicopter32plots_P.Integrator_UpperSat;
  } else if (helicopter32plots_X.Integrator_CSTATE <=
             (helicopter32plots_P.Integrator_LowerSat) ) {
    helicopter32plots_X.Integrator_CSTATE =
      (helicopter32plots_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter32plots_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter32plots_P.elevation_ref_Value -
    helicopter32plots_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter32plots_B.Sum2 = ((helicopter32plots_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter32plots_P.K_ed * helicopter32plots_B.Gain1[5]) +
    helicopter32plots_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter32plots_B.Sum2 - helicopter32plots_B.Sum_k) *
    helicopter32plots_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter32plots_B.K_ei = helicopter32plots_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter32plots_DW.TimeStampA >= helicopter32plots_M->Timing.t[0]) &&
      (helicopter32plots_DW.TimeStampB >= helicopter32plots_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter32plots_DW.TimeStampA;
    lastU = &helicopter32plots_DW.LastUAtTimeA;
    if (helicopter32plots_DW.TimeStampA < helicopter32plots_DW.TimeStampB) {
      if (helicopter32plots_DW.TimeStampB < helicopter32plots_M->Timing.t[0]) {
        rtb_Derivative = helicopter32plots_DW.TimeStampB;
        lastU = &helicopter32plots_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter32plots_DW.TimeStampA >= helicopter32plots_M->Timing.t[0]) {
        rtb_Derivative = helicopter32plots_DW.TimeStampB;
        lastU = &helicopter32plots_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter32plots_B.PitchCounttorad - *lastU) /
      (helicopter32plots_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helicopter32plots_B.Gain_l = helicopter32plots_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter32plots_P.BackmotorSaturation_UpperSat) {
    helicopter32plots_B.BackmotorSaturation =
      helicopter32plots_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter32plots_P.BackmotorSaturation_LowerSat) {
    helicopter32plots_B.BackmotorSaturation =
      helicopter32plots_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter32plots_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helicopter32plots_B.Sum_k + helicopter32plots_B.Sum2) *
    helicopter32plots_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter32plots_P.FrontmotorSaturation_UpperSat) {
    helicopter32plots_B.FrontmotorSaturation =
      helicopter32plots_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter32plots_P.FrontmotorSaturation_LowerSat)
  {
    helicopter32plots_B.FrontmotorSaturation =
      helicopter32plots_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter32plots_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter32plots/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter32plots_DW.HILWriteAnalog_Buffer[0] =
        helicopter32plots_B.FrontmotorSaturation;
      helicopter32plots_DW.HILWriteAnalog_Buffer[1] =
        helicopter32plots_B.BackmotorSaturation;
      result = hil_write_analog(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILWriteAnalog_channels, 2,
        &helicopter32plots_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter32plots_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter32plots_DW.TimeStampA == (rtInf)) {
    helicopter32plots_DW.TimeStampA = helicopter32plots_M->Timing.t[0];
    lastU = &helicopter32plots_DW.LastUAtTimeA;
  } else if (helicopter32plots_DW.TimeStampB == (rtInf)) {
    helicopter32plots_DW.TimeStampB = helicopter32plots_M->Timing.t[0];
    lastU = &helicopter32plots_DW.LastUAtTimeB;
  } else if (helicopter32plots_DW.TimeStampA < helicopter32plots_DW.TimeStampB)
  {
    helicopter32plots_DW.TimeStampA = helicopter32plots_M->Timing.t[0];
    lastU = &helicopter32plots_DW.LastUAtTimeA;
  } else {
    helicopter32plots_DW.TimeStampB = helicopter32plots_M->Timing.t[0];
    lastU = &helicopter32plots_DW.LastUAtTimeB;
  }

  *lastU = helicopter32plots_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter32plots_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter32plots_M->solverInfo);
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
  if (!(++helicopter32plots_M->Timing.clockTick0)) {
    ++helicopter32plots_M->Timing.clockTickH0;
  }

  helicopter32plots_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter32plots_M->solverInfo);

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
    if (!(++helicopter32plots_M->Timing.clockTick1)) {
      ++helicopter32plots_M->Timing.clockTickH1;
    }

    helicopter32plots_M->Timing.t[1] = helicopter32plots_M->Timing.clockTick1 *
      helicopter32plots_M->Timing.stepSize1 +
      helicopter32plots_M->Timing.clockTickH1 *
      helicopter32plots_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter32plots_derivatives(void)
{
  XDot_helicopter32plots_T *_rtXdot;
  _rtXdot = ((XDot_helicopter32plots_T *) helicopter32plots_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter32plots_P.TravelTransferFcn_A *
    helicopter32plots_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter32plots_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter32plots_P.PitchTransferFcn_A *
    helicopter32plots_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter32plots_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter32plots_P.ElevationTransferFcn_A *
    helicopter32plots_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter32plots_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter32plots_X.Integrator_CSTATE <=
            (helicopter32plots_P.Integrator_LowerSat) );
    usat = ( helicopter32plots_X.Integrator_CSTATE >=
            helicopter32plots_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter32plots_B.K_ei > 0)) ||
        (usat && (helicopter32plots_B.K_ei < 0)) ) {
      ((XDot_helicopter32plots_T *) helicopter32plots_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter32plots_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter32plots_T *) helicopter32plots_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter32plots_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter32plots/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter32plots_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter32plots_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter32plots_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
      return;
    }

    if ((helicopter32plots_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter32plots_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter32plots_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter32plots_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter32plots_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter32plots_DW.HILInitialize_Card,
         helicopter32plots_P.HILInitialize_analog_input_chan, 8U,
         &helicopter32plots_DW.HILInitialize_AIMinimums[0],
         &helicopter32plots_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter32plots_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter32plots_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter32plots_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter32plots_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter32plots_DW.HILInitialize_Card,
         helicopter32plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter32plots_DW.HILInitialize_AOMinimums[0],
         &helicopter32plots_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter32plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter32plots_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_analog_output_cha, 8U,
        &helicopter32plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter32plots_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter32plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter32plots_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter32plots_DW.HILInitialize_Card,
         helicopter32plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter32plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter32plots_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter32plots_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter32plots_DW.HILInitialize_Card,
         helicopter32plots_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter32plots_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter32plots_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter32plots_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_encoder_channels, 8U,
        &helicopter32plots_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter32plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter32plots_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter32plots_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter32plots_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter32plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter32plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter32plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter32plots_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter32plots_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter32plots_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter32plots_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter32plots_DW.HILInitialize_Card,
          &helicopter32plots_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter32plots_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter32plots_DW.HILInitialize_Card,
          &helicopter32plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter32plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter32plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter32plots_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter32plots_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter32plots_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter32plots_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter32plots_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *)
        &helicopter32plots_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter32plots_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter32plots_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter32plots_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter32plots_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter32plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter32plots_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter32plots_DW.HILInitialize_POSortedFreqs[0],
        &helicopter32plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter32plots_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter32plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter32plots_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter32plots_DW.HILInitialize_Card,
        helicopter32plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter32plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter32plots_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter32plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter32plots_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter32plots_DW.HILInitialize_Card,
         helicopter32plots_P.HILInitialize_pwm_channels, 8U,
         &helicopter32plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter32plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter32plots_DW.HILInitialize_Card,
       helicopter32plots_P.HILReadEncoderTimebase_samples_,
       helicopter32plots_P.HILReadEncoderTimebase_channels, 3,
       &helicopter32plots_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
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
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625279, 3.1262155534580049,
      3.1033093000299821, 3.0666274151912156, 3.014453922394225,
      2.9456562771176764, 2.8595077632937147, 2.7555515879654058,
      2.6335051104906526, 2.4931956060326308, 2.3345185760651095,
      2.1583782719260478, 1.9678000778704727, 1.7674110445466207,
      1.5625810512761342, 1.3586839679386824, 1.1605920903248914,
      0.9723541967977124, 0.79695012325984016, 0.6364329337495902,
      0.49215959796927389, 0.36485729820996016, 0.25463993267205981,
      0.16109623593017897, 0.083400591252192915, 0.020423624533390444,
      -0.029167512815941796, -0.066823004653065, -0.094039185304793951,
      -0.1123016902197959, -0.12304167596519093, -0.12760358353335186,
      -0.12722285379205736, -0.12301203948381058, -0.115953847476759,
      -0.10689976340781566, -0.0965730445362397, -0.085575006872611575,
      -0.0743936735427865, -0.063413988602610077, -0.052928930978585746,
      -0.043150984592951393, -0.034223531491125968, -0.026231833998759252,
      -0.019213359191285466, -0.013167274288011302, -0.0080630053551753473,
      -0.0038478045493674519, -0.00045331387359379615, 0.0021988530036483609,
      0.0041924638884134014, 0.0056128908698021623, 0.0065441252851624904,
      0.007066500218929463, 0.0072550329904089838, 0.0071783005100552856,
      0.0068977633952542186, 0.0064674596661119139, 0.0059339951101777082,
      0.0053367645211473734, 0.0047083455725129356, 0.0040750147523033907,
      0.0034573422992259205, 0.0028708302460222221, 0.0023265643498006213,
      0.0018318567755375489, 0.0013908618413420552, 0.0010051519083085056,
      0.00067424460561774221, 0.0003960760448359, 0.00016741753288223579,
      -1.5764411858470461E-5, -0.00015800314738642138, -0.00026407916771299959,
      -0.00033881565345682865, -0.0003869178172458332, -0.00041285121080772658,
      -0.00042075405089002664, -0.00041437869183475541, -0.00039705758179327445,
      -0.00037168935026994556, -0.00034074105498252188, -0.0003062630387661211,
      -0.00026991328973278821, -0.00023298864146415992, -0.00019646057949526015,
      -0.00016101382346239441, -0.00012708622092749204, -9.4908810282185824E-5,
      -6.4545177863867975E-5, -3.5929439305405009E-5, -8.9023060705418551E-6,
      1.6755258848342874E-5, 4.1291353364613383E-5, 6.4956154500282379E-5,
      8.7979679297890811E-5, 0.00011055569399486009, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909061995,
      -0.046506351618091045, -0.091625013712091127, -0.14672753935506738,
      -0.20869397118796099, -0.27519058110619471, -0.34459405529584691,
      -0.41582470131323507, -0.48818590989901217, -0.56123801783208871,
      -0.63470811987008446, -0.70456121655624582, -0.76231277622230043,
      -0.801556133295409, -0.81931997308194582, -0.81558833334980685,
      -0.79236751045516318, -0.75295157410871649, -0.70161629415148907,
      -0.64206875804099994, -0.57709334312126559, -0.5092091990372547,
      -0.44086946215160139, -0.37417478696752343, -0.3107825787119442,
      -0.25190786687520988, -0.19836454939732895, -0.15062196734849281,
      -0.10886472260691585, -0.0730500196600078, -0.042959942981580117,
      -0.018247630272643769, 0.0015229189651780454, 0.016843257232987116,
      0.028232768028206272, 0.036216336275773323, 0.041306875486303907,
      0.043992150654512466, 0.044725333319300287, 0.043918739760705677,
      0.041940230496097367, 0.039111785542537421, 0.035709812407301693,
      0.03196678996946687, 0.028073899229895141, 0.024184339613096649,
      0.020417075731343826, 0.016860803223231583, 0.013577962703094623,
      0.010608667508968628, 0.0079744435390601623, 0.0056817079255550461,
      0.0037249376614413112, 0.0020894997350678918, 0.00075413108591808429,
      -0.00030692992141479326, -0.0011221484592042674, -0.00172121491656922,
      -0.0021338582237368218, -0.0023889223561213385, -0.0025136757945377515,
      -0.0025333232808381788, -0.0024706898123098816, -0.0023460482128147926,
      -0.0021770635848864031, -0.0019788302970522894, -0.001763979736781975,
      -0.0015428397321341982, -0.0013236292107630537, -0.0011126742431273688,
      -0.00091463404781465685, -0.000732727778962825, -0.00056895494211180371,
      -0.00042430408130631283, -0.00029894594297531606, -0.00019240865515601823,
      -0.00010373357424757347, -3.1611360329200382E-5, 2.5501436221084938E-5,
      6.9284440165923874E-5, 0.00010147292609331559, 0.00012379318114969466,
      0.00013791206486560333, 0.00014539899613333135, 0.00014769859307451328,
      0.00014611224787559912, 0.00014178702413146296, 0.00013571041013960944,
      0.00012870964258122482, 0.00012145452967327141, 0.00011446295423385185,
      0.00010810853293945261, 0.00010263025967553892, 9.8144378065082008E-5,
      9.4659204542676012E-5, 9.20940991904337E-5, 9.0304058787877108E-5,
      8.9110958714119329E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205861008, 0.22266037932307312, 0.31888147181624243,
      0.38944360631121494, 0.43795507377648224, 0.46997264230352059,
      0.49051724877498, 0.50343100141409236, 0.5114213858593829,
      0.51630439857559984, 0.51925862126752032, 0.49369500885904033,
      0.40816596705876146, 0.27735705984392084, 0.12554803518857205,
      -0.026373804426941864, -0.16411590764817893, -0.27857678423584442,
      -0.3628181526028984, -0.42085924264294045, -0.45922141703380032,
      -0.47977920385573469, -0.48299915977627905, -0.4713724919567861,
      -0.44803191699491934, -0.41610397764324641, -0.37842371849810486,
      -0.33742633592109811, -0.29512425777261325, -0.253124641962517,
      -0.21266516986452458, -0.17465718802102478, -0.1397306911813434,
      -0.10827829968414156, -0.080496712382112937, -0.056424811263810659,
      -0.035977986944966414, -0.018978499319992532, -0.0051818550552452242,
      0.0057006952153475839, 0.013983347843646941, 0.019990368682992504,
      0.024043846827354096, 0.02645425304348089, 0.027513464961974855,
      0.027489921858446689, 0.026625582310723324, 0.025134375864333883,
      0.023201863003346496, 0.020985844389337642, 0.018617688948044882,
      0.016204179667314442, 0.013829704890780045, 0.011558650651908637,
      0.0094378756039130117, 0.0074991739560053856, 0.0057616532742171466,
      0.0042339729232725156, 0.002916405297005134, 0.0018026958727821762,
      0.00088170965649384848, 0.00013886092934043144, -0.00044266943438280026,
      -0.00088091922171180763, -0.0011943188110467246, -0.0014010371685183715,
      -0.0015184817035754469, -0.00156293309435058, -0.0015492962434657841,
      -0.0014909491426511202, -0.0013996724643236184, -0.0012856440340190626,
      -0.0011574838622158488, -0.001022337038652467, -0.0008859834445411937,
      -0.000752964861244103, -0.00062672160478955167, -0.00050973226282778118,
      -0.00040365143331496246, -0.00030944155013681607, -0.000227495925005011,
      -0.00015775103811508151, -9.97868777746858E-5, -5.29147707674965E-5,
      -1.6252672910833967E-5, 1.121168200387841E-5, 3.056902952022862E-5,
      4.2947186894513477E-5, 4.9478751347060973E-5, 5.1276367137421222E-5,
      4.9413784960899296E-5, 4.4910622807850555E-5, 3.8718343149644826E-5,
      3.1704497960415808E-5, 2.46318753877148E-5, 1.812918495638472E-5,
      1.2651321907297626E-5, 8.43237564871279E-6, 5.4485883366234313E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823444032,
      0.46652650905785209, 0.38488436997267739, 0.28224853797988991,
      0.1940458698610692, 0.12807027410815358, 0.082178425885837486,
      0.051655010556449575, 0.03196153778116162, 0.01953205086486786,
      0.011816890767682152, -0.10225444963392009, -0.34211616720111565,
      -0.52323562885936248, -0.60723609862139516, -0.60768735846205568,
      -0.55096841288494824, -0.45784350635066196, -0.3369654734682162,
      -0.23216436016016803, -0.15344869756343943, -0.082231147287737588,
      -0.012879823682177514, 0.046506671277971735, 0.09336229984746712,
      0.12771175740669161, 0.15072103658056607, 0.16398953030802702,
      0.1692083125939394, 0.16799846324038503, 0.16183788839196958,
      0.15203192737399915, 0.1397059873587255, 0.12580956598880741,
      0.11112634920811451, 0.09628760447320911, 0.08178729727537698,
      0.067997950499895529, 0.055186577058989231, 0.043530201082371232,
      0.033130610513197427, 0.024028083357382254, 0.016213912577446359,
      0.0096416248645071843, 0.00423684767397585, -9.41724141126613E-5,
      -0.0034573581908934684, -0.0059648257855577681, -0.0077300514439495385,
      -0.0088640744560354173, -0.00947262176517103, -0.0096540371229217447,
      -0.0094978991061375941, -0.00908421695548563, -0.0084831001919825014,
      -0.0077548065916305071, -0.0069500827271529557, -0.006110721403778524,
      -0.0052702705050695258, -0.0044548376968918313, -0.0036839448651533107,
      -0.0029713949086136683, -0.0023261214548929267, -0.0017529991493160295,
      -0.0012535983573396668, -0.00082687342988658759, -0.00046977814022830245,
      -0.00017780556310053119, 5.4547403539183411E-5, 0.00023338840325865532,
      0.00036510671331000749, 0.00045611372121822345, 0.0005126406872128556,
      0.0005405872942535275, 0.00054541437644509269, 0.00053207433318836264,
      0.00050497302581820535, 0.00046795736784708181, 0.00042432331805127487,
      0.00037683953271258548, 0.0003277825005272203, 0.00027897954755971805,
      0.0002318566413615828, 0.00018748842802875719, 0.00014664839142665014,
      0.00010985741965884952, 7.742939006540083E-5, 4.9512629497139416E-5,
      2.6126257810190002E-5, 7.19046316144098E-6, -7.4503287060876969E-6,
      -1.8012648612194969E-5, -2.4769118632822925E-5, -2.805538075691608E-5,
      -2.8290490290804035E-5, -2.60107617253203E-5, -2.191145219634838E-5,
      -1.6875785034339347E-5, -1.1935149248357435E-5, -8.0237392002931046E-6,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter32plots_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter32plots_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter32plots_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "32.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter32plots_M, "Error creating .mat file 32.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,11,0,"x_state")) {
      rtmSetErrorStatus(helicopter32plots_M,
                        "Error writing mat file header to file 32.mat");
      return;
    }

    helicopter32plots_DW.ToFile_IWORK.Count = 0;
    helicopter32plots_DW.ToFile_IWORK.Decimation = -1;
    helicopter32plots_DW.ToFile_PWORK.FilePtr = fp;
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
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559807444,
      0.52359877559803447, 0.52359877559798307, 0.52359877559791523,
      0.523598775597823, 0.5235987755976923, 0.52359877559749723,
      0.52359877559718215, 0.5235987755966055, 0.52359877559526735,
      0.52359877558925361, 0.38860546189573769, 0.10951698887526899,
      -0.11003833946102244, -0.27691123755291608, -0.39790682731646931,
      -0.47963675127901145, -0.523591909441963, -0.52359877534342569,
      -0.52359877536343258, -0.5235987224145312, -0.503434893916266,
      -0.46497044043306429, -0.4207734586135255, -0.37347862555446432,
      -0.32522986550994049, -0.27772344682220712, -0.23225530276673972,
      -0.18977018817324473, -0.15091074088283757, -0.11606494657763662,
      -0.085410890684983448, -0.05895801601103625, -0.036584388943327439,
      -0.0180697126514252, -0.0031240162861197718, 0.0085879009858575127,
      0.017426078578606414, 0.023758780967378858, 0.027949332008985622,
      0.030346034333594649, 0.031274796827770436, 0.031034093821064389,
      0.029891891654304288, 0.028084198910718453, 0.025814923189764541,
      0.023256747740391895, 0.020552773751695549, 0.017818707168261894,
      0.015145401405518497, 0.01260159841128787, 0.010236739517233149,
      0.0080837440173715247, 0.006161677142781272, 0.0044782499569109881,
      0.00303211167359721, 0.0018149100878760943, 0.00081310836210016517,
      9.5565276842184743E-6, -0.00061517602288165539, -0.001081694820851813,
      -0.0014108848952404617, -0.0016232063484833857, -0.0017381600630382559,
      -0.0017739020793310274, -0.0017469853873323272, -0.001672208684182087,
      -0.0015625529081204714, -0.001429187925983797, -0.0012815335105175714,
      -0.0011273605978063498, -0.000972920686000228, -0.00082309306297666371,
      -0.00068154128644839368, -0.00055087195231493544, -0.00043279025403042583,
      -0.00032824814515925642, -0.00023758206462032373, -0.00016063817112712354,
      -9.6883866645734384E-5, -4.5505078677972753E-5, -5.489330992376534E-6,
      2.4304922656237084E-5, 4.5091880270843033E-5, 5.8113799538615935E-5,
      6.4606933573472E-5, 6.5776522273753516E-5, 6.2779201128075108E-5,
      5.67109984837793E-5, 4.85988379814542E-5, 3.9393103171466169E-5,
      2.995832301200616E-5, 2.10584099342284E-5, 1.3332266149043414E-5,
      7.2554391803702257E-6, 3.0850995699961835E-6, 7.9190498794152718E-7,
      -4.0593744088170127E-17, -1.533214493789752E-17, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter32plots_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter32plots_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter32plots_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter32plots_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter32plots_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter32plots_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter32plots_X.Integrator_CSTATE = helicopter32plots_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter32plots_DW.TimeStampA = (rtInf);
  helicopter32plots_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter32plots_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter32plots/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter32plots_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter32plots_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter32plots_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter32plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter32plots_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter32plots_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter32plots_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter32plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter32plots_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter32plots_DW.HILInitialize_Card
                         , helicopter32plots_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter32plots_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter32plots_DW.HILInitialize_AOVoltages[0]
                         , &helicopter32plots_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter32plots_DW.HILInitialize_Card,
             helicopter32plots_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter32plots_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter32plots_DW.HILInitialize_Card,
            helicopter32plots_P.HILInitialize_pwm_channels,
            num_final_pwm_outputs, &helicopter32plots_DW.HILInitialize_POValues
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter32plots_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter32plots_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter32plots_DW.HILInitialize_Card);
    hil_close(helicopter32plots_DW.HILInitialize_Card);
    helicopter32plots_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter32plots_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "32.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter32plots_M, "Error closing MAT-file 32.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter32plots_M, "Error reopening MAT-file 32.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 11, helicopter32plots_DW.ToFile_IWORK.Count,
           "x_state")) {
        rtmSetErrorStatus(helicopter32plots_M,
                          "Error writing header for x_state to MAT-file 32.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter32plots_M, "Error closing MAT-file 32.mat");
        return;
      }

      helicopter32plots_DW.ToFile_PWORK.FilePtr = (NULL);
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
  helicopter32plots_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter32plots_update();
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
  helicopter32plots_initialize();
}

void MdlTerminate(void)
{
  helicopter32plots_terminate();
}

/* Registration function */
RT_MODEL_helicopter32plots_T *helicopter32plots(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter32plots_P.Integrator_UpperSat = rtInf;
  helicopter32plots_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter32plots_M, 0,
                sizeof(RT_MODEL_helicopter32plots_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter32plots_M->solverInfo,
                          &helicopter32plots_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter32plots_M->solverInfo, &rtmGetTPtr
                (helicopter32plots_M));
    rtsiSetStepSizePtr(&helicopter32plots_M->solverInfo,
                       &helicopter32plots_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter32plots_M->solverInfo,
                 &helicopter32plots_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter32plots_M->solverInfo, (real_T **)
                         &helicopter32plots_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter32plots_M->solverInfo,
      &helicopter32plots_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter32plots_M->solverInfo, (&rtmGetErrorStatus
      (helicopter32plots_M)));
    rtsiSetRTModelPtr(&helicopter32plots_M->solverInfo, helicopter32plots_M);
  }

  rtsiSetSimTimeStep(&helicopter32plots_M->solverInfo, MAJOR_TIME_STEP);
  helicopter32plots_M->ModelData.intgData.f[0] =
    helicopter32plots_M->ModelData.odeF[0];
  helicopter32plots_M->ModelData.contStates = ((real_T *) &helicopter32plots_X);
  rtsiSetSolverData(&helicopter32plots_M->solverInfo, (void *)
                    &helicopter32plots_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter32plots_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter32plots_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter32plots_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter32plots_M->Timing.sampleTimes =
      (&helicopter32plots_M->Timing.sampleTimesArray[0]);
    helicopter32plots_M->Timing.offsetTimes =
      (&helicopter32plots_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter32plots_M->Timing.sampleTimes[0] = (0.0);
    helicopter32plots_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter32plots_M->Timing.offsetTimes[0] = (0.0);
    helicopter32plots_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter32plots_M, &helicopter32plots_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter32plots_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter32plots_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter32plots_M, -1);
  helicopter32plots_M->Timing.stepSize0 = 0.002;
  helicopter32plots_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter32plots_M->Sizes.checksums[0] = (1552757371U);
  helicopter32plots_M->Sizes.checksums[1] = (2338904514U);
  helicopter32plots_M->Sizes.checksums[2] = (546617563U);
  helicopter32plots_M->Sizes.checksums[3] = (3764842860U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter32plots_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter32plots_M->extModeInfo,
      &helicopter32plots_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter32plots_M->extModeInfo,
                        helicopter32plots_M->Sizes.checksums);
    rteiSetTPtr(helicopter32plots_M->extModeInfo, rtmGetTPtr(helicopter32plots_M));
  }

  helicopter32plots_M->solverInfoPtr = (&helicopter32plots_M->solverInfo);
  helicopter32plots_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter32plots_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter32plots_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter32plots_M->ModelData.blockIO = ((void *) &helicopter32plots_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter32plots_B.Gain1[i] = 0.0;
    }

    helicopter32plots_B.FromWorkspace1[0] = 0.0;
    helicopter32plots_B.FromWorkspace1[1] = 0.0;
    helicopter32plots_B.FromWorkspace1[2] = 0.0;
    helicopter32plots_B.FromWorkspace1[3] = 0.0;
    helicopter32plots_B.TravelCounttorad = 0.0;
    helicopter32plots_B.Gain = 0.0;
    helicopter32plots_B.Sum1 = 0.0;
    helicopter32plots_B.Gain_d = 0.0;
    helicopter32plots_B.PitchCounttorad = 0.0;
    helicopter32plots_B.Gain_i = 0.0;
    helicopter32plots_B.Gain_b = 0.0;
    helicopter32plots_B.ElevationCounttorad = 0.0;
    helicopter32plots_B.Gain_e = 0.0;
    helicopter32plots_B.Sum = 0.0;
    helicopter32plots_B.Gain_dg = 0.0;
    helicopter32plots_B.FromWorkspace = 0.0;
    helicopter32plots_B.Sum_g = 0.0;
    helicopter32plots_B.TmpSignalConversionAtToWorkspac[0] = 0.0;
    helicopter32plots_B.TmpSignalConversionAtToWorkspac[1] = 0.0;
    helicopter32plots_B.Sum_k = 0.0;
    helicopter32plots_B.Sum2 = 0.0;
    helicopter32plots_B.K_ei = 0.0;
    helicopter32plots_B.Gain_l = 0.0;
    helicopter32plots_B.BackmotorSaturation = 0.0;
    helicopter32plots_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter32plots_M->ModelData.defaultParam = ((real_T *)&helicopter32plots_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter32plots_X;
    helicopter32plots_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter32plots_X, 0,
                  sizeof(X_helicopter32plots_T));
  }

  /* states (dwork) */
  helicopter32plots_M->ModelData.dwork = ((void *) &helicopter32plots_DW);
  (void) memset((void *)&helicopter32plots_DW, 0,
                sizeof(DW_helicopter32plots_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter32plots_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter32plots_DW.TimeStampA = 0.0;
  helicopter32plots_DW.LastUAtTimeA = 0.0;
  helicopter32plots_DW.TimeStampB = 0.0;
  helicopter32plots_DW.LastUAtTimeB = 0.0;
  helicopter32plots_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter32plots_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter32plots_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter32plots_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter32plots_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter32plots_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter32plots_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter32plots_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter32plots_M->Sizes.numBlocks = (64);/* Number of blocks */
  helicopter32plots_M->Sizes.numBlockIO = (22);/* Number of block outputs */
  helicopter32plots_M->Sizes.numBlockPrms = (146);/* Sum of parameter "widths" */
  return helicopter32plots_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
