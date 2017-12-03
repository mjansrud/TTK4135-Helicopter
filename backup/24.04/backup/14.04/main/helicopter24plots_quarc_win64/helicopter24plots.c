/*
 * helicopter24plots.c
 *
 * Code generation for model "helicopter24plots".
 *
 * Model version              : 1.187
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Sun Apr 23 12:53:46 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter24plots.h"
#include "helicopter24plots_private.h"
#include "helicopter24plots_dt.h"

/* Block signals (auto storage) */
B_helicopter24plots_T helicopter24plots_B;

/* Continuous states */
X_helicopter24plots_T helicopter24plots_X;

/* Block states (auto storage) */
DW_helicopter24plots_T helicopter24plots_DW;

/* Real-time model */
RT_MODEL_helicopter24plots_T helicopter24plots_M_;
RT_MODEL_helicopter24plots_T *const helicopter24plots_M = &helicopter24plots_M_;

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
  helicopter24plots_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter24plots_output(void)
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
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* set solver stop time */
    if (!(helicopter24plots_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter24plots_M->solverInfo,
                            ((helicopter24plots_M->Timing.clockTickH0 + 1) *
        helicopter24plots_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter24plots_M->solverInfo,
                            ((helicopter24plots_M->Timing.clockTick0 + 1) *
        helicopter24plots_M->Timing.stepSize0 +
        helicopter24plots_M->Timing.clockTickH0 *
        helicopter24plots_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter24plots_M)) {
    helicopter24plots_M->Timing.t[0] = rtsiGetT(&helicopter24plots_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter24plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter24plots_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter24plots_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter24plots_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter24plots_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter24plots_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter24plots_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter24plots_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter24plots_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter24plots_M->Timing.t[0];

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

    helicopter24plots_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter24plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter24plots_B.FromWorkspace1[0])[elIdx] =
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
            (&helicopter24plots_B.FromWorkspace1[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter24plots_B.TravelCounttorad =
      helicopter24plots_P.TravelCounttorad_Gain * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter24plots_B.Gain = helicopter24plots_P.Gain_Gain *
      helicopter24plots_B.TravelCounttorad;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24plots_P.TravelTransferFcn_C *
    helicopter24plots_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter24plots_P.TravelTransferFcn_D *
    helicopter24plots_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter24plots_B.Gain_d = helicopter24plots_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter24plots_B.PitchCounttorad =
      helicopter24plots_P.PitchCounttorad_Gain * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter24plots_B.Gain_i = helicopter24plots_P.Gain_Gain_a *
      helicopter24plots_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24plots_P.PitchTransferFcn_C *
    helicopter24plots_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter24plots_P.PitchTransferFcn_D *
    helicopter24plots_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter24plots_B.Gain_b = helicopter24plots_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter24plots_B.ElevationCounttorad =
      helicopter24plots_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter24plots_B.Gain_e = helicopter24plots_P.Gain_Gain_lv *
      helicopter24plots_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter24plots_B.Sum = helicopter24plots_B.Gain_e +
      helicopter24plots_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24plots_P.ElevationTransferFcn_C *
    helicopter24plots_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter24plots_P.ElevationTransferFcn_D *
    helicopter24plots_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter24plots_B.Gain_dg = helicopter24plots_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter24plots_B.Gain1[0] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Gain;
  helicopter24plots_B.Gain1[1] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Gain_d;
  helicopter24plots_B.Gain1[2] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Gain_i;
  helicopter24plots_B.Gain1[3] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Gain_b;
  helicopter24plots_B.Gain1[4] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Sum;
  helicopter24plots_B.Gain1[5] = helicopter24plots_P.Gain1_Gain *
    helicopter24plots_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    rtb_TmpSignalConversionAtToFile[0] = helicopter24plots_B.FromWorkspace1[0];
    rtb_TmpSignalConversionAtToFile[1] = helicopter24plots_B.FromWorkspace1[1];
    rtb_TmpSignalConversionAtToFile[2] = helicopter24plots_B.FromWorkspace1[2];
    rtb_TmpSignalConversionAtToFile[3] = helicopter24plots_B.FromWorkspace1[3];
    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i + 4] = helicopter24plots_B.Gain1[i];
    }

    /* End of SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter24plots_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter24plots_DW.ToFile_IWORK.Count*11)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter24plots_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[11];
          helicopter24plots_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter24plots_M->Timing.t[1];
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
            rtmSetErrorStatus(helicopter24plots_M,
                              "Error writing to MAT-file 24.mat");
            return;
          }

          if (((++helicopter24plots_DW.ToFile_IWORK.Count)*11)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file 24.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter24plots_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter24plots_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter24plots_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter24plots_M->Timing.t[0];

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

    helicopter24plots_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Backgain = pDataValues[currTimeIndex];
        } else {
          rtb_Backgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Backgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter24plots_B.Sum_k = ((rtb_Backgain - helicopter24plots_B.Gain1[2]) *
    helicopter24plots_P.K_pp - helicopter24plots_P.K_pd *
    helicopter24plots_B.Gain1[3]) + helicopter24plots_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter24plots_X.Integrator_CSTATE >=
      helicopter24plots_P.Integrator_UpperSat ) {
    helicopter24plots_X.Integrator_CSTATE =
      helicopter24plots_P.Integrator_UpperSat;
  } else if (helicopter24plots_X.Integrator_CSTATE <=
             (helicopter24plots_P.Integrator_LowerSat) ) {
    helicopter24plots_X.Integrator_CSTATE =
      (helicopter24plots_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter24plots_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter24plots_P.elevation_ref_Value -
    helicopter24plots_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter24plots_B.Sum2 = ((helicopter24plots_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter24plots_P.K_ed * helicopter24plots_B.Gain1[5]) +
    helicopter24plots_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter24plots_B.Sum2 - helicopter24plots_B.Sum_k) *
    helicopter24plots_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter24plots_B.K_ei = helicopter24plots_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter24plots_DW.TimeStampA >= helicopter24plots_M->Timing.t[0]) &&
      (helicopter24plots_DW.TimeStampB >= helicopter24plots_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter24plots_DW.TimeStampA;
    lastU = &helicopter24plots_DW.LastUAtTimeA;
    if (helicopter24plots_DW.TimeStampA < helicopter24plots_DW.TimeStampB) {
      if (helicopter24plots_DW.TimeStampB < helicopter24plots_M->Timing.t[0]) {
        rtb_Derivative = helicopter24plots_DW.TimeStampB;
        lastU = &helicopter24plots_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter24plots_DW.TimeStampA >= helicopter24plots_M->Timing.t[0]) {
        rtb_Derivative = helicopter24plots_DW.TimeStampB;
        lastU = &helicopter24plots_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter24plots_B.PitchCounttorad - *lastU) /
      (helicopter24plots_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter24plots_B.Gain_l = helicopter24plots_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter24plots_P.BackmotorSaturation_UpperSat) {
    helicopter24plots_B.BackmotorSaturation =
      helicopter24plots_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter24plots_P.BackmotorSaturation_LowerSat) {
    helicopter24plots_B.BackmotorSaturation =
      helicopter24plots_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter24plots_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helicopter24plots_B.Sum_k + helicopter24plots_B.Sum2) *
    helicopter24plots_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter24plots_P.FrontmotorSaturation_UpperSat) {
    helicopter24plots_B.FrontmotorSaturation =
      helicopter24plots_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter24plots_P.FrontmotorSaturation_LowerSat)
  {
    helicopter24plots_B.FrontmotorSaturation =
      helicopter24plots_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter24plots_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter24plots/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter24plots_DW.HILWriteAnalog_Buffer[0] =
        helicopter24plots_B.FrontmotorSaturation;
      helicopter24plots_DW.HILWriteAnalog_Buffer[1] =
        helicopter24plots_B.BackmotorSaturation;
      result = hil_write_analog(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILWriteAnalog_channels, 2,
        &helicopter24plots_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter24plots_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter24plots_DW.TimeStampA == (rtInf)) {
    helicopter24plots_DW.TimeStampA = helicopter24plots_M->Timing.t[0];
    lastU = &helicopter24plots_DW.LastUAtTimeA;
  } else if (helicopter24plots_DW.TimeStampB == (rtInf)) {
    helicopter24plots_DW.TimeStampB = helicopter24plots_M->Timing.t[0];
    lastU = &helicopter24plots_DW.LastUAtTimeB;
  } else if (helicopter24plots_DW.TimeStampA < helicopter24plots_DW.TimeStampB)
  {
    helicopter24plots_DW.TimeStampA = helicopter24plots_M->Timing.t[0];
    lastU = &helicopter24plots_DW.LastUAtTimeA;
  } else {
    helicopter24plots_DW.TimeStampB = helicopter24plots_M->Timing.t[0];
    lastU = &helicopter24plots_DW.LastUAtTimeB;
  }

  *lastU = helicopter24plots_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter24plots_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter24plots_M->solverInfo);
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
  if (!(++helicopter24plots_M->Timing.clockTick0)) {
    ++helicopter24plots_M->Timing.clockTickH0;
  }

  helicopter24plots_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter24plots_M->solverInfo);

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
    if (!(++helicopter24plots_M->Timing.clockTick1)) {
      ++helicopter24plots_M->Timing.clockTickH1;
    }

    helicopter24plots_M->Timing.t[1] = helicopter24plots_M->Timing.clockTick1 *
      helicopter24plots_M->Timing.stepSize1 +
      helicopter24plots_M->Timing.clockTickH1 *
      helicopter24plots_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter24plots_derivatives(void)
{
  XDot_helicopter24plots_T *_rtXdot;
  _rtXdot = ((XDot_helicopter24plots_T *) helicopter24plots_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter24plots_P.TravelTransferFcn_A *
    helicopter24plots_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter24plots_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter24plots_P.PitchTransferFcn_A *
    helicopter24plots_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter24plots_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter24plots_P.ElevationTransferFcn_A *
    helicopter24plots_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter24plots_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter24plots_X.Integrator_CSTATE <=
            (helicopter24plots_P.Integrator_LowerSat) );
    usat = ( helicopter24plots_X.Integrator_CSTATE >=
            helicopter24plots_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter24plots_B.K_ei > 0)) ||
        (usat && (helicopter24plots_B.K_ei < 0)) ) {
      ((XDot_helicopter24plots_T *) helicopter24plots_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter24plots_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter24plots_T *) helicopter24plots_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter24plots_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter24plots/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter24plots_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter24plots_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter24plots_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
      return;
    }

    if ((helicopter24plots_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter24plots_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter24plots_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter24plots_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter24plots_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter24plots_DW.HILInitialize_Card,
         helicopter24plots_P.HILInitialize_analog_input_chan, 8U,
         &helicopter24plots_DW.HILInitialize_AIMinimums[0],
         &helicopter24plots_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter24plots_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter24plots_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter24plots_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter24plots_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter24plots_DW.HILInitialize_Card,
         helicopter24plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter24plots_DW.HILInitialize_AOMinimums[0],
         &helicopter24plots_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter24plots_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_analog_output_cha, 8U,
        &helicopter24plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter24plots_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter24plots_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter24plots_DW.HILInitialize_Card,
         helicopter24plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter24plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter24plots_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter24plots_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter24plots_DW.HILInitialize_Card,
         helicopter24plots_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter24plots_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter24plots_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter24plots_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_encoder_channels, 8U,
        &helicopter24plots_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter24plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter24plots_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter24plots_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter24plots_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter24plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter24plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter24plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter24plots_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter24plots_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter24plots_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter24plots_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter24plots_DW.HILInitialize_Card,
          &helicopter24plots_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter24plots_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter24plots_DW.HILInitialize_Card,
          &helicopter24plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter24plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter24plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter24plots_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter24plots_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter24plots_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter24plots_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter24plots_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *)
        &helicopter24plots_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter24plots_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter24plots_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter24plots_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter24plots_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter24plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24plots_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter24plots_DW.HILInitialize_POSortedFreqs[0],
        &helicopter24plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24plots_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24plots_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter24plots_DW.HILInitialize_Card,
        helicopter24plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter24plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter24plots_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24plots_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter24plots_DW.HILInitialize_Card,
         helicopter24plots_P.HILInitialize_pwm_channels, 8U,
         &helicopter24plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter24plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter24plots_DW.HILInitialize_Card,
       helicopter24plots_P.HILReadEncoderTimebase_samples_,
       helicopter24plots_P.HILReadEncoderTimebase_channels, 3,
       &helicopter24plots_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
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
      3.1415926535897931, 3.1378421413510189, 3.12621555344675,
      3.1033093000193528, 3.066627415181765, 3.0144539223866809,
      2.945656277112954, 2.8595077632929788, 2.7555515879701851,
      2.6335051105028726, 2.4931956060551728, 2.3345185761018206,
      2.1574113215284938, 1.9618364841188816, 1.7516021208794057,
      1.5345792352435073, 1.3182739937579158, 1.1085622010868816,
      0.90963662698130432, 0.72433188725680941, 0.55449774646178518,
      0.40131147431604658, 0.26550851373440115, 0.14754172782854258,
      0.047686794539707915, -0.033890031482728, -0.097089262730253728,
      -0.14185175494675403, -0.16814274063400478, -0.1773935551104629,
      -0.17248741351119173, -0.15710545488210803, -0.13505558539896353,
      -0.10978494359791881, -0.084102278247095041, -0.060073535490778575,
      -0.039038964138158154, -0.021704012966765392, -0.0082663574444696141,
      0.0014478425902237821, 0.0078559378157679466, 0.011509048217845154,
      0.013005057027954235, 0.012924963188563121, 0.011790752295337263,
      0.010041559374534449, 0.0080243077816887949, 0.0059949768546746609,
      0.0041269789773387177, 0.0025236582303639165, 0.001232540141534777,
      0.00025958155609148638, -0.00041776350379775723, -0.00083936346922836513,
      -0.0010534309932294286, -0.0011097771703820595, -0.001055132710569317,
      -0.00093025225104187762, -0.00076846184611268015, -0.00059530015714915718,
      -0.00042892658963453585, -0.00028101219168505333, -0.00015788152306709028,
      -6.1728352023554528E-5, 8.2203544526594072E-6, 5.4672427418278445E-5,
      8.1413564872156823E-5, 9.2639896583660849E-5, 9.2473950518694875E-5,
      8.464460812551682E-5, 7.2304280676270126E-5, 5.7953816238311565E-5,
      4.3446273198789706E-5, 3.004352717375679E-5, 1.8503758403434684E-5,
      9.1824483999964016E-6, 2.1340517778924414E-6, -2.7943623002618614E-6,
      -5.8826960810305788E-6, -7.4730356465872176E-6, -7.9208032807171376E-6,
      -7.5607317212340552E-6, -6.6856708992200183E-6, -5.535816535127497E-6,
      -4.2958670725820952E-6, -3.0977672802775478E-6, -2.0269939965227003E-6,
      -1.1307095043680898E-6, -4.2649601991255812E-7, 8.9247696707428208E-8,
      4.33843641947182E-7, 6.3213577001080327E-7, 7.1167922058413627E-7,
      6.9929276609479541E-7, 6.1885012979501724E-7, 4.9012573177316218E-7,
      3.2847899639220526E-7, 1.4514845134319154E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048903438231,
      -0.04650635161167551, -0.0916250137041888, -0.14672753934495231,
      -0.20869397117493577, -0.27519058108950983, -0.34459405527450071,
      -0.4158247012857742, -0.48818590986384847, -0.5612380177853995,
      -0.63470811980800956, -0.70842901828790694, -0.78229934963304859,
      -0.84093745295250388, -0.86809154253819421, -0.86522096593696651,
      -0.83884717067873615, -0.79570229641690926, -0.74121895889257994,
      -0.67933656317469693, -0.61274508857755439, -0.54321184232118158,
      -0.47186714361803439, -0.39941973314993867, -0.32630730408434372,
      -0.25279692498470296, -0.17904996886060129, -0.10516394274360295,
      -0.037003257900432546, 0.019624566402484611, 0.061527834521734721,
      0.088199477937978013, 0.10108256720957888, 0.10273066140869502,
      0.09611497103066588, 0.084138285415881628, 0.069339804690971021,
      0.053750622094583075, 0.038856800144173559, 0.025632380907576621,
      0.0146124416137088, 0.0059840352458362953, -0.00032037535216448591,
      -0.0045368435675034648, -0.0069967716778112877, -0.0080690063659826471,
      -0.0081173237026565658, -0.0074719915039438058, -0.0064132829824992375,
      -0.0051644723499165894, -0.0038918343363731943, -0.0027093802341570062,
      -0.0016863998563224629, -0.00085627009060428537, -0.00022538470321055581,
      0.00021857784465093778, 0.00049952184350972621, 0.00064716162511675843,
      0.00069264676125406, 0.000665494275458454, 0.00059165759719789833,
      0.00049252267987182064, 0.00038461268957411135, 0.00027979483130482412,
      0.00018580829726244451, 0.00010696455521548189, 4.4905332245984418E-5,
      -6.6377885989551984E-7, -3.1317364172743846E-5, -4.9361304397018439E-5,
      -5.7401852351865881E-5, -5.80301667581191E-5, -5.36109787001633E-5,
      -4.6159069681320077E-5, -3.7285234613784773E-5, -2.8193581088447487E-5,
      -1.971365091264886E-5, -1.2353329723106521E-5, -6.3613528622582045E-6,
      -1.7910651365513257E-6, 1.4402916379006786E-6, 3.5002486880245004E-6,
      4.5994228563384382E-6, 4.9598032501499613E-6, 4.7924045691865387E-6,
      4.28309853498774E-6, 3.5851433685867939E-6, 2.8168593377904777E-6,
      2.0629802664482966E-6, 1.3783891809273663E-6, 7.9317391222283606E-7,
      3.1817920226168323E-7, -4.9540417989012308E-8, -3.2176514523076186E-7,
      -5.14892192119069E-7, -6.4658154155547641E-7, -7.3331678022770376E-7,
      -7.8888955582973282E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205731112, 0.22266037931775992, 0.31888147180601706,
      0.38944360629585911, 0.43795507375619691, 0.46997264227793856,
      0.49051724874231778, 0.50343100137115937, 0.51142138580522456,
      0.51630439849442578, 0.51925862115906185, 0.52103115473418249,
      0.5220872891538485, 0.41443171900968473, 0.19191473441494833,
      -0.020288139077372289, -0.1863999121204071, -0.304931493253539,
      -0.38506742117064863, -0.43736113861064763, -0.47064310962274314,
      -0.49143442817931388, -0.50423708223738672, -0.51203062784082853,
      -0.51673072530191733, -0.51954328415619933, -0.52121532021083339,
      -0.52219821380656917, -0.48173368832385738, -0.40022383469243783,
      -0.29615629523600756, -0.18850498914611905, -0.09105275462061356,
      -0.011648100353647219, 0.046757172906502983, 0.084646639751391237,
      0.10459000988072469, 0.11017838872093146, 0.10526384524952433,
      0.093465144451936238, 0.0778847221587777, 0.060982280819769,
      0.044557166306234568, 0.029800386979233043, 0.017385832380304053,
      0.0075781452649664963, 0.00034148848238179291, -0.0045609615128757441,
      -0.0074825474839330308, -0.0088261166010489209, -0.0089945194294497131,
      -0.0083571339875177057, -0.0072300346103045844, -0.00586704013799048,
      -0.00445885696804535, -0.0031377577284019561, -0.0019856048846179719,
      -0.0010434615890567997, -0.00032147157010661316, 0.00019190339972161545,
      0.00052184944228734923, 0.00070064773401784693, 0.0007626665983773229,
      0.00074081259015741857, 0.0006642609272557268, 0.0005572374568158109,
      0.00043861088637365143, 0.00032206507367801948, 0.00021664783395832616,
      0.0001275276783883861, 5.6827522355737518E-5, 4.4406865382860309E-6,
      -3.1233133120411813E-5, -5.2667246605446915E-5, -6.2716876851045034E-5,
      -6.425633452086537E-5, -5.9932907517669941E-5, -5.2019938807465806E-5,
      -4.2349003694711587E-5, -3.2301047863606627E-5, -2.2837995268940845E-5,
      -1.4558989402129367E-5, -7.7685429387771459E-6, -2.5470307030821953E-6,
      1.1831102263499675E-6, 3.5995808166645529E-6, 4.9328808335876478E-6,
      5.4299383869086474E-6, 5.32812963507514E-6, 4.838428613980331E-6,
      4.1360782946527341E-6, 3.3570814812438829E-6, 2.59890211801673E-6,
      1.9239806859453483E-6, 1.3649484904470119E-6, 9.30730309339465E-7,
      6.130118234764176E-7, 3.9276715342242723E-7, 2.4660785250980917E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4241150082335135,
      0.46652650904719517, 0.38488436995842851, 0.28224853796476823,
      0.19404586984675123, 0.12807027409236652, 0.082178425862916821,
      0.051655010520766285, 0.031961537741660793, 0.019532050762204708,
      0.011816890663944193, 0.0070901343058826539, 0.0042245376840640367,
      -0.43062228057125523, -0.89006793837354559, -0.84881149396388245,
      -0.66444709216673936, -0.47412632452712772, -0.3205437116630383,
      -0.20917486975459612, -0.13312788404298195, -0.083165274220882879,
      -0.051210616226891534, -0.03117418240836738, -0.018800389838955395,
      -0.01125023541172789, -0.0066881442131362491, -0.0039315743775431532,
      0.16185810193624711, 0.32603941453107815, 0.4162701578311212,
      0.43060522436495396, 0.38980893810742195, 0.31761861707326533,
      0.23362109304600079, 0.15155786738495294, 0.079773480522733753,
      0.022353515366227081, -0.019658173880228547, -0.047194803184952394,
      -0.062321689167234143, -0.067609765350634857, -0.065700458048737753,
      -0.059027117302606132, -0.049658218390316, -0.03923074845595026,
      -0.028946627124938847, -0.019609799975630178, -0.011686343878829178,
      -0.0053742764630635882, -0.00067361130820320152, 0.0025495417731279963,
      0.0045083975142524556, 0.0054519778946563848, 0.00563273268518049,
      0.0052843969639735427, 0.004608611380535906, 0.0037685731876446575,
      0.0028879600812007143, 0.0020534998847128827, 0.0013197841756629033,
      0.00071519317232195915, 0.00024807546283787203, -8.7416027479648851E-5,
      -0.00030620664620679885, -0.00042809387635969509, -0.00047450627636866935,
      -0.00046618324538255928, -0.000421668953478805, -0.00035648061687979188,
      -0.00028280061873062597, -0.00020954733786983758, -0.00014269527323482303,
      -8.573644854017208E-5, -4.0198515582424109E-5, -6.1578252793130249E-6,
      1.7293713412750077E-5, 3.165188024078491E-5, 3.8683745850985223E-5,
      4.0191828724388187E-5, 3.785221577863147E-5, 3.3116028867214262E-5,
      2.7161791253377238E-5, 2.0886054342748152E-5, 1.4920569117697003E-5,
      9.6658877612266932E-6, 5.3332054676607311E-6, 1.9882356132523487E-6,
      -4.0722960736567649E-7, -1.9587986844108843E-6, -2.8093958773420356E-6,
      -3.1159818536670537E-6, -3.0327120529402615E-6, -2.6996803283171757E-6,
      -2.2361233820249942E-6, -1.7368673244618372E-6, -1.270868543483838E-6,
      -8.8097328024761025E-7, -5.8463180368212108E-7, -3.7660457464037824E-7,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter24plots_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter24plots_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter24plots_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "24.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter24plots_M, "Error creating .mat file 24.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,11,0,"x_state")) {
      rtmSetErrorStatus(helicopter24plots_M,
                        "Error writing mat file header to file 24.mat");
      return;
    }

    helicopter24plots_DW.ToFile_IWORK.Count = 0;
    helicopter24plots_DW.ToFile_IWORK.Decimation = -1;
    helicopter24plots_DW.ToFile_PWORK.FilePtr = fp;
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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559054624,
      0.52359877557927514, 0.52359877557440859, 0.52359877556969248,
      0.52359877556578793, 0.523598775558614, 0.52359877554486389,
      0.52359877552437639, 0.52359877551358225, 0.52359877542446642,
      0.52359877541369826, 0.52359877532730015, 0.52359877455241466,
      -0.011122838283745327, -0.52359877297118829, -0.52359877497968288,
      -0.52359877518407749, -0.52359877526911769, -0.5235987753490311,
      -0.52359877546730726, -0.5235987755809971, -0.52359877542724587,
      -0.52359877481985828, -0.52359877483429929, -0.52359877557439738,
      -0.52359877158489321, -0.52359877552021794, -0.52357138287359373,
      -0.32090512345880695, -0.1396628998157132, -0.008071692801445882,
      0.079996175896428615, 0.13192816854925973, 0.15549220006087247,
      0.1581561729519966, 0.14661826764078845, 0.12653210349832827,
      0.1023949944539761, 0.077560966150216987, 0.054340134644648125,
      0.034150007158661518, 0.017690334033699927, 0.005119930160587309,
      -0.0037795358500980632, -0.0094619777245251255, -0.012501991566172661,
      -0.013507442044638895, -0.013057826670831359, -0.011665133599255403,
      -0.0097531122289753, -0.007650675467583419, -0.0055953722811867678,
      -0.00374335482194393, -0.0021828893722446604, -0.00094912733764423963,
      -3.8492371823169595E-5, 0.00057839201871832859, 0.00094583738855001274,
      0.001114521976587308, 0.0011351851661689282, 0.0010543730994386257,
      0.00091192160946938989, 0.00073981998641784851, 0.00056210010608771408,
      0.00039542579798373592, 0.00025010491534385493, 0.00013130179936353991,
      4.02836017470902E-5, -2.4414587454151713E-5, -6.5976563842823512E-5,
      -8.847878567453589E-5, -9.6259088101424364E-5, -9.3469445925510339E-5,
      -8.3790086395863712E-5, -7.027618940236047E-5, -5.5306726335128553E-5,
      -4.0606394035353849E-5, -2.7314971509089596E-5, -1.6082842126391908E-5,
      -7.1761742650977283E-6, -5.7983326727809967E-7, 3.909819333158307E-6,
      6.6066757269656854E-6, 7.8729493107581766E-6, 8.0733815346187421E-6,
      7.544112493813453E-6, 6.5739890307531092E-6, 5.3957662568062319E-6,
      4.1846610451155158E-6, 3.0619376663582372E-6, 2.1015606611585195E-6,
      1.3383704338514843E-6, 7.7667245817733731E-7, 3.9854569976898791E-7,
      1.7154383509723319E-7, 5.5758653561578452E-8, 1.0397044756678316E-8,
      3.26782992949579E-13, 2.9800866950168836E-13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    helicopter24plots_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter24plots_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter24plots_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter24plots_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter24plots_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter24plots_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter24plots_X.Integrator_CSTATE = helicopter24plots_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter24plots_DW.TimeStampA = (rtInf);
  helicopter24plots_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter24plots_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter24plots/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter24plots_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter24plots_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter24plots_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24plots_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter24plots_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter24plots_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter24plots_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24plots_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter24plots_DW.HILInitialize_Card
                         , helicopter24plots_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter24plots_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter24plots_DW.HILInitialize_AOVoltages[0]
                         , &helicopter24plots_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter24plots_DW.HILInitialize_Card,
             helicopter24plots_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter24plots_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter24plots_DW.HILInitialize_Card,
            helicopter24plots_P.HILInitialize_pwm_channels,
            num_final_pwm_outputs, &helicopter24plots_DW.HILInitialize_POValues
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24plots_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter24plots_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter24plots_DW.HILInitialize_Card);
    hil_close(helicopter24plots_DW.HILInitialize_Card);
    helicopter24plots_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter24plots_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "24.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter24plots_M, "Error closing MAT-file 24.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter24plots_M, "Error reopening MAT-file 24.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 11, helicopter24plots_DW.ToFile_IWORK.Count,
           "x_state")) {
        rtmSetErrorStatus(helicopter24plots_M,
                          "Error writing header for x_state to MAT-file 24.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter24plots_M, "Error closing MAT-file 24.mat");
        return;
      }

      helicopter24plots_DW.ToFile_PWORK.FilePtr = (NULL);
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
  helicopter24plots_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter24plots_update();
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
  helicopter24plots_initialize();
}

void MdlTerminate(void)
{
  helicopter24plots_terminate();
}

/* Registration function */
RT_MODEL_helicopter24plots_T *helicopter24plots(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter24plots_P.Integrator_UpperSat = rtInf;
  helicopter24plots_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter24plots_M, 0,
                sizeof(RT_MODEL_helicopter24plots_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter24plots_M->solverInfo,
                          &helicopter24plots_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter24plots_M->solverInfo, &rtmGetTPtr
                (helicopter24plots_M));
    rtsiSetStepSizePtr(&helicopter24plots_M->solverInfo,
                       &helicopter24plots_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter24plots_M->solverInfo,
                 &helicopter24plots_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter24plots_M->solverInfo, (real_T **)
                         &helicopter24plots_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter24plots_M->solverInfo,
      &helicopter24plots_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter24plots_M->solverInfo, (&rtmGetErrorStatus
      (helicopter24plots_M)));
    rtsiSetRTModelPtr(&helicopter24plots_M->solverInfo, helicopter24plots_M);
  }

  rtsiSetSimTimeStep(&helicopter24plots_M->solverInfo, MAJOR_TIME_STEP);
  helicopter24plots_M->ModelData.intgData.f[0] =
    helicopter24plots_M->ModelData.odeF[0];
  helicopter24plots_M->ModelData.contStates = ((real_T *) &helicopter24plots_X);
  rtsiSetSolverData(&helicopter24plots_M->solverInfo, (void *)
                    &helicopter24plots_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter24plots_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter24plots_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter24plots_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter24plots_M->Timing.sampleTimes =
      (&helicopter24plots_M->Timing.sampleTimesArray[0]);
    helicopter24plots_M->Timing.offsetTimes =
      (&helicopter24plots_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter24plots_M->Timing.sampleTimes[0] = (0.0);
    helicopter24plots_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter24plots_M->Timing.offsetTimes[0] = (0.0);
    helicopter24plots_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter24plots_M, &helicopter24plots_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter24plots_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter24plots_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter24plots_M, -1);
  helicopter24plots_M->Timing.stepSize0 = 0.002;
  helicopter24plots_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter24plots_M->Sizes.checksums[0] = (2738611720U);
  helicopter24plots_M->Sizes.checksums[1] = (2145499042U);
  helicopter24plots_M->Sizes.checksums[2] = (3489053730U);
  helicopter24plots_M->Sizes.checksums[3] = (1170973073U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter24plots_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter24plots_M->extModeInfo,
      &helicopter24plots_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter24plots_M->extModeInfo,
                        helicopter24plots_M->Sizes.checksums);
    rteiSetTPtr(helicopter24plots_M->extModeInfo, rtmGetTPtr(helicopter24plots_M));
  }

  helicopter24plots_M->solverInfoPtr = (&helicopter24plots_M->solverInfo);
  helicopter24plots_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter24plots_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter24plots_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter24plots_M->ModelData.blockIO = ((void *) &helicopter24plots_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter24plots_B.Gain1[i] = 0.0;
    }

    helicopter24plots_B.FromWorkspace1[0] = 0.0;
    helicopter24plots_B.FromWorkspace1[1] = 0.0;
    helicopter24plots_B.FromWorkspace1[2] = 0.0;
    helicopter24plots_B.FromWorkspace1[3] = 0.0;
    helicopter24plots_B.TravelCounttorad = 0.0;
    helicopter24plots_B.Gain = 0.0;
    helicopter24plots_B.Gain_d = 0.0;
    helicopter24plots_B.PitchCounttorad = 0.0;
    helicopter24plots_B.Gain_i = 0.0;
    helicopter24plots_B.Gain_b = 0.0;
    helicopter24plots_B.ElevationCounttorad = 0.0;
    helicopter24plots_B.Gain_e = 0.0;
    helicopter24plots_B.Sum = 0.0;
    helicopter24plots_B.Gain_dg = 0.0;
    helicopter24plots_B.Sum_k = 0.0;
    helicopter24plots_B.Sum2 = 0.0;
    helicopter24plots_B.K_ei = 0.0;
    helicopter24plots_B.Gain_l = 0.0;
    helicopter24plots_B.BackmotorSaturation = 0.0;
    helicopter24plots_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter24plots_M->ModelData.defaultParam = ((real_T *)&helicopter24plots_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter24plots_X;
    helicopter24plots_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter24plots_X, 0,
                  sizeof(X_helicopter24plots_T));
  }

  /* states (dwork) */
  helicopter24plots_M->ModelData.dwork = ((void *) &helicopter24plots_DW);
  (void) memset((void *)&helicopter24plots_DW, 0,
                sizeof(DW_helicopter24plots_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24plots_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter24plots_DW.TimeStampA = 0.0;
  helicopter24plots_DW.LastUAtTimeA = 0.0;
  helicopter24plots_DW.TimeStampB = 0.0;
  helicopter24plots_DW.LastUAtTimeB = 0.0;
  helicopter24plots_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter24plots_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter24plots_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter24plots_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter24plots_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter24plots_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter24plots_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter24plots_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter24plots_M->Sizes.numBlocks = (56);/* Number of blocks */
  helicopter24plots_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helicopter24plots_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helicopter24plots_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
