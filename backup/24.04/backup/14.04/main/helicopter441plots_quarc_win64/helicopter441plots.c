/*
 * helicopter441plots.c
 *
 * Code generation for model "helicopter441plots".
 *
 * Model version              : 1.179
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Apr 05 11:07:56 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter441plots.h"
#include "helicopter441plots_private.h"
#include "helicopter441plots_dt.h"

/* Block signals (auto storage) */
B_helicopter441plots_T helicopter441plots_B;

/* Continuous states */
X_helicopter441plots_T helicopter441plots_X;

/* Block states (auto storage) */
DW_helicopter441plots_T helicopter441plots_DW;

/* Real-time model */
RT_MODEL_helicopter441plots_T helicopter441plots_M_;
RT_MODEL_helicopter441plots_T *const helicopter441plots_M =
  &helicopter441plots_M_;

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
  helicopter441plots_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter441plots_output(void)
{
  /* local block i/o variables */
  real_T rtb_Derivative;
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[12];
  real_T *lastU;
  real_T lastTime;
  int32_T i;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* set solver stop time */
    if (!(helicopter441plots_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter441plots_M->solverInfo,
                            ((helicopter441plots_M->Timing.clockTickH0 + 1) *
        helicopter441plots_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter441plots_M->solverInfo,
                            ((helicopter441plots_M->Timing.clockTick0 + 1) *
        helicopter441plots_M->Timing.stepSize0 +
        helicopter441plots_M->Timing.clockTickH0 *
        helicopter441plots_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter441plots_M)) {
    helicopter441plots_M->Timing.t[0] = rtsiGetT
      (&helicopter441plots_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter441plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter441plots_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter441plots_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter441plots_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter441plots_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter441plots_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter441plots_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter441plots_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter441plots_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter441plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[160]) {
      currTimeIndex = 159;
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

    helicopter441plots_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&helicopter441plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 161;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&helicopter441plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 161;
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
            (&helicopter441plots_B.FromWorkspace1[0])[elIdx] = (real_T)
              rtInterpolate(d1, d2, f1, f2);
            pDataValues += 161;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter441plots_B.TravelCounttorad =
      helicopter441plots_P.TravelCounttorad_Gain * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter441plots_B.Gain = helicopter441plots_P.Gain_Gain *
      helicopter441plots_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    helicopter441plots_B.Sum1 = helicopter441plots_P.travel_offsetdeg_Value +
      helicopter441plots_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter441plots_P.TravelTransferFcn_C *
    helicopter441plots_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter441plots_P.TravelTransferFcn_D *
    helicopter441plots_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter441plots_B.Gain_d = helicopter441plots_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter441plots_B.PitchCounttorad =
      helicopter441plots_P.PitchCounttorad_Gain * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter441plots_B.Gain_i = helicopter441plots_P.Gain_Gain_a *
      helicopter441plots_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter441plots_P.PitchTransferFcn_C *
    helicopter441plots_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter441plots_P.PitchTransferFcn_D *
    helicopter441plots_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter441plots_B.Gain_b = helicopter441plots_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter441plots_B.ElevationCounttorad =
      helicopter441plots_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter441plots_B.Gain_e = helicopter441plots_P.Gain_Gain_lv *
      helicopter441plots_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter441plots_B.Sum = helicopter441plots_B.Gain_e +
      helicopter441plots_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter441plots_P.ElevationTransferFcn_C *
    helicopter441plots_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter441plots_P.ElevationTransferFcn_D *
    helicopter441plots_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter441plots_B.Gain_dg = helicopter441plots_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter441plots_B.Gain1[0] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Sum1;
  helicopter441plots_B.Gain1[1] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Gain_d;
  helicopter441plots_B.Gain1[2] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Gain_i;
  helicopter441plots_B.Gain1[3] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Gain_b;
  helicopter441plots_B.Gain1[4] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Sum;
  helicopter441plots_B.Gain1[5] = helicopter441plots_P.Gain1_Gain *
    helicopter441plots_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i] = helicopter441plots_B.FromWorkspace1[i];
    }

    for (i = 0; i < 6; i++) {
      rtb_TmpSignalConversionAtToFile[i + 6] = helicopter441plots_B.Gain1[i];
    }

    /* End of SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter441plots_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter441plots_DW.ToFile_IWORK.Count*13)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter441plots_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[13];
          helicopter441plots_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter441plots_M->Timing.t[1];
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
            rtmSetErrorStatus(helicopter441plots_M,
                              "Error writing to MAT-file 441.mat");
            return;
          }

          if (((++helicopter441plots_DW.ToFile_IWORK.Count)*13)+1 >= 100000000)
          {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file 441.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter441plots_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter441plots_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter441plots_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter441plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[160]) {
      currTimeIndex = 159;
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

    helicopter441plots_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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
        pDataValues += 161;
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
  helicopter441plots_B.Sum_k = ((rtb_Backgain - helicopter441plots_B.Gain1[2]) *
    helicopter441plots_P.K_pp - helicopter441plots_P.K_pd *
    helicopter441plots_B.Gain1[3]) + helicopter441plots_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter441plots_X.Integrator_CSTATE >=
      helicopter441plots_P.Integrator_UpperSat ) {
    helicopter441plots_X.Integrator_CSTATE =
      helicopter441plots_P.Integrator_UpperSat;
  } else if (helicopter441plots_X.Integrator_CSTATE <=
             (helicopter441plots_P.Integrator_LowerSat) ) {
    helicopter441plots_X.Integrator_CSTATE =
      (helicopter441plots_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter441plots_X.Integrator_CSTATE;

  /* FromWorkspace: '<Root>/From Workspace2' */
  {
    real_T *pDataValues = (real_T *)
      helicopter441plots_DW.FromWorkspace2_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter441plots_DW.FromWorkspace2_PWORK.TimePtr;
    int_T currTimeIndex = helicopter441plots_DW.FromWorkspace2_IWORK.PrevIndex;
    real_T t = helicopter441plots_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[160]) {
      currTimeIndex = 159;
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

    helicopter441plots_DW.FromWorkspace2_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Derivative = pDataValues[currTimeIndex];
        } else {
          rtb_Derivative = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Derivative = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 161;
      }
    }
  }

  /* Sum: '<S3>/Sum' */
  rtb_Derivative -= helicopter441plots_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter441plots_B.Sum2 = ((helicopter441plots_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter441plots_P.K_ed * helicopter441plots_B.Gain1[5]) +
    helicopter441plots_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter441plots_B.Sum2 - helicopter441plots_B.Sum_k) *
    helicopter441plots_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter441plots_B.K_ei = helicopter441plots_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter441plots_DW.TimeStampA >= helicopter441plots_M->Timing.t[0]) &&
      (helicopter441plots_DW.TimeStampB >= helicopter441plots_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    lastTime = helicopter441plots_DW.TimeStampA;
    lastU = &helicopter441plots_DW.LastUAtTimeA;
    if (helicopter441plots_DW.TimeStampA < helicopter441plots_DW.TimeStampB) {
      if (helicopter441plots_DW.TimeStampB < helicopter441plots_M->Timing.t[0])
      {
        lastTime = helicopter441plots_DW.TimeStampB;
        lastU = &helicopter441plots_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter441plots_DW.TimeStampA >= helicopter441plots_M->Timing.t[0])
      {
        lastTime = helicopter441plots_DW.TimeStampB;
        lastU = &helicopter441plots_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter441plots_B.PitchCounttorad - *lastU) /
      (helicopter441plots_M->Timing.t[0] - lastTime);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter441plots_B.Gain_l = helicopter441plots_P.Gain_Gain_a1 *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter441plots_P.BackmotorSaturation_UpperSat) {
    helicopter441plots_B.BackmotorSaturation =
      helicopter441plots_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter441plots_P.BackmotorSaturation_LowerSat) {
    helicopter441plots_B.BackmotorSaturation =
      helicopter441plots_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter441plots_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  lastTime = (helicopter441plots_B.Sum_k + helicopter441plots_B.Sum2) *
    helicopter441plots_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (lastTime > helicopter441plots_P.FrontmotorSaturation_UpperSat) {
    helicopter441plots_B.FrontmotorSaturation =
      helicopter441plots_P.FrontmotorSaturation_UpperSat;
  } else if (lastTime < helicopter441plots_P.FrontmotorSaturation_LowerSat) {
    helicopter441plots_B.FrontmotorSaturation =
      helicopter441plots_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter441plots_B.FrontmotorSaturation = lastTime;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter441plots/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter441plots_DW.HILWriteAnalog_Buffer[0] =
        helicopter441plots_B.FrontmotorSaturation;
      helicopter441plots_DW.HILWriteAnalog_Buffer[1] =
        helicopter441plots_B.BackmotorSaturation;
      result = hil_write_analog(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILWriteAnalog_channels, 2,
        &helicopter441plots_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter441plots_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter441plots_DW.TimeStampA == (rtInf)) {
    helicopter441plots_DW.TimeStampA = helicopter441plots_M->Timing.t[0];
    lastU = &helicopter441plots_DW.LastUAtTimeA;
  } else if (helicopter441plots_DW.TimeStampB == (rtInf)) {
    helicopter441plots_DW.TimeStampB = helicopter441plots_M->Timing.t[0];
    lastU = &helicopter441plots_DW.LastUAtTimeB;
  } else if (helicopter441plots_DW.TimeStampA < helicopter441plots_DW.TimeStampB)
  {
    helicopter441plots_DW.TimeStampA = helicopter441plots_M->Timing.t[0];
    lastU = &helicopter441plots_DW.LastUAtTimeA;
  } else {
    helicopter441plots_DW.TimeStampB = helicopter441plots_M->Timing.t[0];
    lastU = &helicopter441plots_DW.LastUAtTimeB;
  }

  *lastU = helicopter441plots_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter441plots_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter441plots_M->solverInfo);
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
  if (!(++helicopter441plots_M->Timing.clockTick0)) {
    ++helicopter441plots_M->Timing.clockTickH0;
  }

  helicopter441plots_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter441plots_M->solverInfo);

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
    if (!(++helicopter441plots_M->Timing.clockTick1)) {
      ++helicopter441plots_M->Timing.clockTickH1;
    }

    helicopter441plots_M->Timing.t[1] = helicopter441plots_M->Timing.clockTick1 *
      helicopter441plots_M->Timing.stepSize1 +
      helicopter441plots_M->Timing.clockTickH1 *
      helicopter441plots_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter441plots_derivatives(void)
{
  XDot_helicopter441plots_T *_rtXdot;
  _rtXdot = ((XDot_helicopter441plots_T *)
             helicopter441plots_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter441plots_P.TravelTransferFcn_A *
    helicopter441plots_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter441plots_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter441plots_P.PitchTransferFcn_A *
    helicopter441plots_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter441plots_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter441plots_P.ElevationTransferFcn_A *
    helicopter441plots_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter441plots_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter441plots_X.Integrator_CSTATE <=
            (helicopter441plots_P.Integrator_LowerSat) );
    usat = ( helicopter441plots_X.Integrator_CSTATE >=
            helicopter441plots_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter441plots_B.K_ei > 0)) ||
        (usat && (helicopter441plots_B.K_ei < 0)) ) {
      ((XDot_helicopter441plots_T *) helicopter441plots_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter441plots_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter441plots_T *) helicopter441plots_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter441plots_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter441plots/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter441plots_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter441plots_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter441plots_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
      return;
    }

    if ((helicopter441plots_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_analog_inpu_m && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter441plots_DW.HILInitialize_AIMinimums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter441plots_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter441plots_DW.HILInitialize_AIMaximums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter441plots_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_analog_input_chan, 8U,
         &helicopter441plots_DW.HILInitialize_AIMinimums[0],
         &helicopter441plots_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_analog_outp_b && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter441plots_DW.HILInitialize_AOMinimums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter441plots_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter441plots_DW.HILInitialize_AOMaximums
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter441plots_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter441plots_DW.HILInitialize_AOMinimums[0],
         &helicopter441plots_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_analog_outp_j && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter441plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter441plots_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILInitialize_analog_output_cha, 8U,
        &helicopter441plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter441plots_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter441plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter441plots_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_analog_output_cha, 8U,
         &helicopter441plots_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_encoder_par_m && is_switching))
    {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter441plots_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter441plots_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter441plots_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_encoder_cou_k && is_switching))
    {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter441plots_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter441plots_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILInitialize_encoder_channels, 8U,
        &helicopter441plots_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_pwm_params__f && is_switching))
    {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter441plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter441plots_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter441plots_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter441plots_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter441plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter441plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter441plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter441plots_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter441plots_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter441plots_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter441plots_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter441plots_DW.HILInitialize_Card,
          &helicopter441plots_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter441plots_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter441plots_DW.HILInitialize_Card,
          &helicopter441plots_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter441plots_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter441plots_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter441plots_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter441plots_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            helicopter441plots_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter441plots_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            helicopter441plots_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_pwm_channels, 8U,
         (t_pwm_configuration *)
         &helicopter441plots_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *) &helicopter441plots_DW.HILInitialize_POAlignValues
         [0],
         (t_pwm_polarity *) &helicopter441plots_DW.HILInitialize_POPolarityVals
         [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter441plots_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter441plots_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter441plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter441plots_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter441plots_DW.HILInitialize_POSortedFreqs[0],
        &helicopter441plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter441plots_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_pwm_outputs_g && is_switching))
    {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter441plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter441plots_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter441plots_DW.HILInitialize_Card,
        helicopter441plots_P.HILInitialize_pwm_channels, 8U,
        &helicopter441plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }

    if (helicopter441plots_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter441plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter441plots_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter441plots_DW.HILInitialize_Card,
         helicopter441plots_P.HILInitialize_pwm_channels, 8U,
         &helicopter441plots_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter441plots/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter441plots_DW.HILInitialize_Card,
       helicopter441plots_P.HILReadEncoderTimebase_samples_,
       helicopter441plots_P.HILReadEncoderTimebase_channels, 3,
       &helicopter441plots_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
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
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0,
      36.25, 36.5, 36.75, 37.0, 37.25, 37.5, 37.75, 38.0, 38.25, 38.5, 38.75,
      39.0, 39.25, 39.5, 39.75, 40.0 } ;

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
      2.5292010906318554, 1.935861443396339, 1.5131752902008089,
      1.21491638044382, 1.057920541228291, 0.97049792076076236,
      0.9183221489465424, 0.872245758821696, 0.818846844953982,
      0.75761853401164059, 0.69070483050500386, 0.62045775266950032,
      0.54907275824179591, 0.47846594836847961, 0.41021761332278706,
      0.34556278031215248, 0.28540980043106234, 0.23037426916099621,
      0.18081926148169919, 0.13689627655756131, 0.098584059744561969,
      0.0657239438450172, 0.038051376076873492, 0.015223560893734044,
      -0.003156630400354239, -0.017520132825795322, -0.02831407953107266,
      -0.035987516318884638, -0.040980109030208989, -0.043713001560459937,
      -0.0445816210910268, -0.043951607607498937, -0.042156083290333921,
      -0.039494420547850059, -0.0362324835830042, -0.032603722603231113,
      -0.028811606651270555, -0.025031343514603781, -0.021412775322283941,
      -0.018083334407660662, -0.015151126331109941, -0.012708616442066575,
      -0.01083584236000959, -0.0094416439639205072, -0.0084015337406797026,
      -0.0076061513626637994, -0.0069728579420073079, -0.0064438322014218116,
      -0.0059804689992299679, -0.0055577936159615345, -0.0051600119157288555,
      -0.0047773052016341373, -0.00440365063309681, -0.0040353932414395636,
      -0.00367033447028956, -0.00330716153572716, -0.0029450954386272534,
      -0.0025836763973080749, -0.0022226342989174592, -0.0018618111033409047,
      -0.0015011146757620089, -0.0011404914747276725, -0.00077991047564329534,
      -0.00041935374766349586, -5.8810951778772428E-5, 0.00030172386069196773,
      0.00066225410598457315, 0.0010227817423329967, 0.0013833078903821453,
      0.0017438331905685262, 0.0021043580082249886, 0.0024648825516646912,
      0.0028254069394042522, 0.0031859312388253804, 0.0035464554881472713,
      0.0039069797091506441, 0.0042675039141263065, 0.0046280281100185908,
      0.0049885523007996858, 0.0053490764886885605, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.14372722710909097, -0.17694571407297521, -0.13765228850408257,
      -0.11541761506840911, -0.109473894749103, -0.13248348866401863,
      -0.17327272462014168, -0.21287272080920297, -0.24489105381953996,
      -0.26765448211629284, -0.28098830888122966, -0.28553997770191036,
      -0.28242723949292337, -0.27299334018279015, -0.25861933204230009,
      -0.24061191952429126, -0.22014212508036152, -0.19822003071698743,
      -0.17569193969665628, -0.15324886725186374, -0.13144046359816136,
      -0.11069027107252297, -0.091311260732764007, -0.073520765176146835,
      -0.057454009701953773, -0.04317578682078789, -0.03069374715121314,
      -0.019970370845308408, -0.010931570120970317, -0.0034744781224925175,
      0.0025200539340837855, 0.00718209726879829, 0.010646650969916752,
      0.013047747859436439, 0.014515043919119879, 0.015168463807591194,
      0.015121052546922174, 0.014474272768967855, 0.01331776365866244,
      0.011728832306107098, 0.0097700395560519453, 0.007491096328217736,
      0.005576793584519936, 0.004160440893002537, 0.0031815295120852915,
      0.0025331736824097161, 0.0021161029624012073, 0.0018534528088853749,
      0.0016907015330227854, 0.0015911268007672756, 0.001530826856371641,
      0.0014946182741720174, 0.0014730295667034227, 0.0014602350845099085,
      0.0014526917382903174, 0.0014482643883749484, 0.0014456761652823973,
      0.001444168393649758, 0.0014432927821528515, 0.001442785710533506,
      0.0014424928042082756, 0.0014423239963028112, 0.0014422269118809972,
      0.0014421711833604977, 0.00144213925000464, 0.0014421209811444483,
      0.0014421105453980669, 0.0014421045922448587, 0.001442101200646052,
      0.0014420992707098626, 0.0014420981737659595, 0.0014420975508979102,
      0.0014420971976099478, 0.0014420969973731673, 0.0014420968840005756,
      0.0014420968198713318, 0.0014420967836082192, 0.0014420967631238381,
      0.0014420967515593706, 0.0014420967450357344, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.0036299999179726389, 0.091977675136779458, 0.19689726429817345,
      0.27692652629216452, 0.3188343059003651, 0.31716469860180058,
      0.282123094040267, 0.22637236155615989, 0.16088454568477034,
      0.094238406000879918, 0.032169456881855574, -0.021999644834657609,
      -0.066675197298447855, -0.10158999976954332, -0.12726951422840765,
      -0.14467268923975093, -0.15493699040419462, -0.15921994333136957,
      -0.15861906450828531, -0.15413346788839871, -0.14665443579078677,
      -0.1369634437867919, -0.12573642799035667, -0.11355369142495667,
      -0.10091302613947027, -0.088218289204380435, -0.075788728223071089,
      -0.063882791391591143, -0.052703878209652638, -0.04236706312433506,
      -0.032949541747009779, -0.024486142365990705, -0.016970035780249741,
      -0.010370288155921403, -0.0046181222158667434, 0.0003350846821730427,
      0.0045711924506756929, 0.0081737646945383054, 0.011229959949469263,
      0.013843999048972475, 0.016106700352953182, 0.01352956067593954,
      0.010010239886755144, 0.0069185717709677238, 0.0045823313823185668,
      0.0029476965598419979, 0.001856310972676688, 0.0011502638606240491,
      0.00070375617742358758, 0.00042617697732668546, 0.00025590843138521342,
      0.00015258073888669543, 9.0426512291254832E-5, 5.3313489984345404E-5,
      3.1290819012826627E-5, 1.8292570204155618E-5, 1.0656354475827363E-5,
      6.1884874356398304E-6, 3.5837889566334511E-6, 2.0701503991861243E-6,
      1.1930692865963968E-6, 6.8615571737390615E-7, 3.9386783091901305E-7,
      2.2569250470266326E-7, 1.2911673373367425E-7, 7.3756431127232346E-8,
      4.2074263773912383E-8, 2.3970368600196734E-8, 1.3639938204398771E-8,
      7.7528947079188609E-9, 4.4021021856613152E-9, 2.497065158037681E-9,
      1.4151352664125094E-9, 8.0128727029011289E-10, 4.5333720698113956E-10,
      2.5628172417065822E-10, 1.4477556984500081E-10, 8.1727836230534361E-11,
      4.6106025439528573E-11, 2.5993958252889828E-11, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.39551023768503196, 0.42180010443609295, 0.31484670429984957,
      0.15535120565320007, -0.010152796929453676, -0.14049864539405346,
      -0.22301612637247409, -0.26195150376538184, -0.266584560826183,
      -0.24827579648499687, -0.21667640686603606, -0.17870220985517357,
      -0.13965920988435079, -0.1027180578354497, -0.0696127000453735,
      -0.041057204657740785, -0.017131811708677145, 0.0024035152923249127,
      0.017942386479566988, 0.029916128390449069, 0.038763968015970554,
      0.044908063185733524, 0.048730946261589206, 0.050562661141929907,
      0.050778947740353415, 0.04971824392523972, 0.047623747325918536,
      0.044715652727762956, 0.041347260341244287, 0.037670085509304639,
      0.033853597524066276, 0.030064426342954691, 0.026398990497340145,
      0.023008663760217769, 0.019812827592177347, 0.016944431074011983,
      0.01441028897545437, 0.0122247810197235, 0.010456156398033611,
      0.0090508052159423667, -0.010308558708071998, -0.014077283156701302,
      -0.012366672463180425, -0.00934496155461195, -0.0065385392899329235,
      -0.0043655423486911238, -0.0028241884483472053, -0.0017860307327438952,
      -0.0011103168003708641, -0.00068107418374332843, -0.00041331077003274093,
      -0.00024861690642463203, -0.00014845208917230618, -8.8090683872221647E-5,
      -5.1992995254445595E-5, -3.0544862931294879E-5, -1.7871468157617157E-5,
      -1.0418793916496336E-5, -6.05455421724926E-6, -3.5083244909015116E-6,
      -2.0276542719891591E-6, -1.1691515536066152E-6, -6.7270129386215414E-7,
      -3.8630306736406134E-7, -2.2144123065661961E-7, -1.267286804351929E-7,
      -7.2415584412423827E-8, -4.1321714009105709E-8, -2.3548170721834743E-8,
      -1.3403166066701794E-8, -7.6201631193636947E-9, -4.3277209375611048E-9,
      -2.4553949086798879E-9, -1.391799652752708E-9, -7.8822236542134531E-10,
      -4.4602553852777334E-10, -2.5219089722557756E-10, -1.4248726251049589E-10,
      -8.0448288611960571E-11, -4.5390715031471693E-11, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.036738921175778261, 0.077207194227187567, 0.08744468007135088,
      0.094652317012421827, 0.10073233246412827, 0.10828932785461459,
      0.11564262846493818, 0.12250789724961078, 0.12887470230529174,
      0.13472995197590468, 0.14008494075031838, 0.14495856975531538,
      0.14937516346961449, 0.15336265604781862, 0.15695125847460964,
      0.16017237889595615, 0.1630577982742536, 0.16563905997588022,
      0.16794695115170585, 0.17001118189525954, 0.17186012432036504,
      0.17352064277575305, 0.17501801312942611, 0.17637589392477507,
      0.1776163053710148, 0.17875965315458284, 0.17982441261412993,
      0.18082950294053723, 0.18179735636000219, 0.18274296492006351,
      0.18368155597054062, 0.18462843586973013, 0.18559863560644346,
      0.1866065785155035, 0.18766629788063593, 0.18879193667524335,
      0.1899966443309283, 0.19129132609083652, 0.19268490532433361,
      0.19418364673944785, 0.19580976363751454, 0.19504879471305761,
      0.19247548504283005, 0.18855542684846627, 0.18366539150903574,
      0.17810887414616774, 0.17212867395986867, 0.16591732408769155,
      0.15962586618122346, 0.15337121825256778, 0.14724234608305242,
      0.14130542299376803, 0.13560813538692831, 0.1301832676049802,
      0.12505167938362713, 0.12022477195892002, 0.11570652425868776,
      0.11149516817195952, 0.10758456131602338, 0.1039653067301195,
      0.1006256612809555, 0.097552268067971418, 0.094730742595267675,
      0.0921461377879121, 0.089783308949692164, 0.08762719638271671,
      0.0856630405270619, 0.083876542054890438, 0.08225397730295031,
      0.0807822776947343, 0.079449080341512143, 0.0782427558094108,
      0.077152419962684932, 0.0761679933865151, 0.07528103324178638,
      0.074488814276873955, 0.073796871191141727, 0.073203522612115185,
      0.07264234062102809, 0.071832515681412273, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.037694070543609688, 0.034029105462999268, 0.030986379853135942,
      0.030310162424603597, 0.029100172134655832, 0.027369801849881931,
      0.025464835594176616, 0.023479061552588047, 0.02143042296485238,
      0.019425482906277179, 0.017495393595197181, 0.015662165229973872,
      0.013940489672462609, 0.012339673283193946, 0.010864655271945056,
      0.0095170360470867334, 0.0082959399427494657, 0.0071983895195965614,
      0.0062201029052594732, 0.0053557393489818564, 0.0045992685459696082,
      0.0039443311520273644, 0.0033844496981355167, 0.0029130625761569448,
      0.0025237053235873238, 0.0022086539562180597, 0.0019696862819838321,
      0.0018208677857198636, 0.0017324478231235129, 0.0017053873483338254,
      0.0017400281203021576, 0.0018352991102544474, 0.0019888039657776587,
      0.0021990201204349913, 0.0024664393539814583, 0.0027872292015600738,
      0.003152709598972295, 0.003554996633535196, 0.0039820927561896981,
      0.0044898768254196314, -0.0050959957515252308, -0.012380712727645255,
      -0.017794494458426638, -0.02169320227183252, -0.024371552756002463,
      -0.026073799691300709, -0.027002122345698012, -0.027323399143214713,
      -0.027174880688856356, -0.026668990135330764, -0.0258973858780283,
      -0.024934396301849262, -0.023839920168136887, -0.022661872615272985,
      -0.021438245084374905, -0.020198837491183882, -0.0189667124091538,
      -0.01775941370318871, -0.016589985783614165, -0.01546782428578823,
      -0.014399384393344003, -0.013388769101748866, -0.012438216368606606,
      -0.01154850123631619, -0.010719266571155088, -0.0099492939804970947,
      -0.0092367246948985085, -0.0085792386899722888, -0.0079741990361423647,
      -0.0074187673699878907, -0.0069099954507343478, -0.0064448969629468328,
      -0.0060205026839051014, -0.0056339119264569547, -0.0052829764251106541,
      -0.0049760380525168844, -0.0047805018174130577, -0.0048950766503647594,
      -0.0056526858130419152, -0.0073769277018600954, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter441plots_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter441plots_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter441plots_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "441.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter441plots_M, "Error creating .mat file 441.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,13,0,"x_state")) {
      rtmSetErrorStatus(helicopter441plots_M,
                        "Error writing mat file header to file 441.mat");
      return;
    }

    helicopter441plots_DW.ToFile_IWORK.Count = 0;
    helicopter441plots_DW.ToFile_IWORK.Decimation = -1;
    helicopter441plots_DW.ToFile_PWORK.FilePtr = fp;
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
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0,
      36.25, 36.5, 36.75, 37.0, 37.25, 37.5, 37.75, 38.0, 38.25, 38.5, 38.75,
      39.0, 39.25, 39.5, 39.75, 40.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.48893682730641569, 0.46928215786367439, 0.42860703735126521,
      0.35150075076547371, 0.24903441979375346, 0.14773940498980162,
      0.05928957126358142, -0.013736617091009306, -0.0704046082004288,
      -0.11271710842886491, -0.14261199796837665, -0.16170038159785372,
      -0.17235642125297146, -0.17624573661575024, -0.17485036305574994,
      -0.16936326071384575, -0.16075428336679132, -0.14985464909161167,
      -0.13736557919515502, -0.12390068211014588, -0.10997006643540216,
      -0.0959980575387661, -0.0823459870670184, -0.069329555748664826,
      -0.057105936330929315, -0.045801483977489552, -0.035561594224880443,
      -0.026463693661131656, -0.018357241924346422, -0.011302199845232471,
      -0.0052231261204636822, -1.2422499419510158E-5, 0.00439354635381797,
      0.0081765873191132037, 0.011249417049208272, 0.013854900394989763,
      0.016033659827917976, 0.017884466675585579, 0.01957336752918307,
      0.021112910043470184, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter441plots_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter441plots_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter441plots_DW.FromWorkspace_IWORK.PrevIndex = 0;
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
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0,
      36.25, 36.5, 36.75, 37.0, 37.25, 37.5, 37.75, 38.0, 38.25, 38.5, 38.75,
      39.0, 39.25, 39.5, 39.75, 40.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.52359877554874712, 0.14444303760782315, 0.14015894787507946,
      0.13653050045201526, 0.13359500450443354, 0.1308813166908607,
      0.1287625229257538, 0.1270864881125584, 0.12500739099943953,
      0.12389608497969525, 0.12294635660970279, 0.12214709825800711,
      0.12148795630613293, 0.12096573355969754, 0.12057565533829091,
      0.12031425804973178, 0.12017976197698944, 0.12016687356259478,
      0.12027443274221956, 0.12049924868790683, 0.12083873637774167,
      0.12129198832577354, 0.12185865065363831, 0.12253685356038886,
      0.12332567442067065, 0.12420123462369947, 0.12530518466691423,
      0.12685982791353417, 0.12823940292184488, 0.12983884585645247,
      0.13166720997774195, 0.1337186148305117, 0.13598301325451168,
      0.13847997327824021, 0.14125020298898905, 0.14424059335725267,
      0.14737151570806312, 0.15062781029260919, 0.1539021283883926,
      0.15819789150479621, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter441plots_DW.FromWorkspace2_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter441plots_DW.FromWorkspace2_PWORK.DataPtr = (void *) pDataValues0;
    helicopter441plots_DW.FromWorkspace2_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter441plots_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter441plots_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter441plots_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter441plots_X.Integrator_CSTATE = helicopter441plots_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter441plots_DW.TimeStampA = (rtInf);
  helicopter441plots_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter441plots_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter441plots/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter441plots_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter441plots_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter441plots_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_analog_outp_c && is_switching))
    {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter441plots_DW.HILInitialize_AOVoltages
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter441plots_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter441plots_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter441plots_P.HILInitialize_set_pwm_outputs_p && is_switching))
    {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter441plots_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter441plots_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter441plots_DW.HILInitialize_Card
                         , helicopter441plots_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter441plots_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter441plots_DW.HILInitialize_AOVoltages[0]
                         , &helicopter441plots_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter441plots_DW.HILInitialize_Card,
             helicopter441plots_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter441plots_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter441plots_DW.HILInitialize_Card,
            helicopter441plots_P.HILInitialize_pwm_channels,
            num_final_pwm_outputs,
            &helicopter441plots_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter441plots_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter441plots_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter441plots_DW.HILInitialize_Card);
    hil_close(helicopter441plots_DW.HILInitialize_Card);
    helicopter441plots_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter441plots_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "441.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter441plots_M, "Error closing MAT-file 441.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter441plots_M,
                          "Error reopening MAT-file 441.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 13,
           helicopter441plots_DW.ToFile_IWORK.Count, "x_state")) {
        rtmSetErrorStatus(helicopter441plots_M,
                          "Error writing header for x_state to MAT-file 441.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter441plots_M, "Error closing MAT-file 441.mat");
        return;
      }

      helicopter441plots_DW.ToFile_PWORK.FilePtr = (NULL);
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
  helicopter441plots_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter441plots_update();
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
  helicopter441plots_initialize();
}

void MdlTerminate(void)
{
  helicopter441plots_terminate();
}

/* Registration function */
RT_MODEL_helicopter441plots_T *helicopter441plots(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter441plots_P.Integrator_UpperSat = rtInf;
  helicopter441plots_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter441plots_M, 0,
                sizeof(RT_MODEL_helicopter441plots_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter441plots_M->solverInfo,
                          &helicopter441plots_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter441plots_M->solverInfo, &rtmGetTPtr
                (helicopter441plots_M));
    rtsiSetStepSizePtr(&helicopter441plots_M->solverInfo,
                       &helicopter441plots_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter441plots_M->solverInfo,
                 &helicopter441plots_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter441plots_M->solverInfo, (real_T **)
                         &helicopter441plots_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter441plots_M->solverInfo,
      &helicopter441plots_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter441plots_M->solverInfo, (&rtmGetErrorStatus
                           (helicopter441plots_M)));
    rtsiSetRTModelPtr(&helicopter441plots_M->solverInfo, helicopter441plots_M);
  }

  rtsiSetSimTimeStep(&helicopter441plots_M->solverInfo, MAJOR_TIME_STEP);
  helicopter441plots_M->ModelData.intgData.f[0] =
    helicopter441plots_M->ModelData.odeF[0];
  helicopter441plots_M->ModelData.contStates = ((real_T *) &helicopter441plots_X);
  rtsiSetSolverData(&helicopter441plots_M->solverInfo, (void *)
                    &helicopter441plots_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter441plots_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter441plots_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter441plots_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter441plots_M->Timing.sampleTimes =
      (&helicopter441plots_M->Timing.sampleTimesArray[0]);
    helicopter441plots_M->Timing.offsetTimes =
      (&helicopter441plots_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter441plots_M->Timing.sampleTimes[0] = (0.0);
    helicopter441plots_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter441plots_M->Timing.offsetTimes[0] = (0.0);
    helicopter441plots_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter441plots_M, &helicopter441plots_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter441plots_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter441plots_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter441plots_M, -1);
  helicopter441plots_M->Timing.stepSize0 = 0.002;
  helicopter441plots_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter441plots_M->Sizes.checksums[0] = (2353948514U);
  helicopter441plots_M->Sizes.checksums[1] = (2790631413U);
  helicopter441plots_M->Sizes.checksums[2] = (3857666905U);
  helicopter441plots_M->Sizes.checksums[3] = (1408898436U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter441plots_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter441plots_M->extModeInfo,
      &helicopter441plots_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter441plots_M->extModeInfo,
                        helicopter441plots_M->Sizes.checksums);
    rteiSetTPtr(helicopter441plots_M->extModeInfo, rtmGetTPtr
                (helicopter441plots_M));
  }

  helicopter441plots_M->solverInfoPtr = (&helicopter441plots_M->solverInfo);
  helicopter441plots_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter441plots_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter441plots_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter441plots_M->ModelData.blockIO = ((void *) &helicopter441plots_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter441plots_B.FromWorkspace1[i] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      helicopter441plots_B.Gain1[i] = 0.0;
    }

    helicopter441plots_B.TravelCounttorad = 0.0;
    helicopter441plots_B.Gain = 0.0;
    helicopter441plots_B.Sum1 = 0.0;
    helicopter441plots_B.Gain_d = 0.0;
    helicopter441plots_B.PitchCounttorad = 0.0;
    helicopter441plots_B.Gain_i = 0.0;
    helicopter441plots_B.Gain_b = 0.0;
    helicopter441plots_B.ElevationCounttorad = 0.0;
    helicopter441plots_B.Gain_e = 0.0;
    helicopter441plots_B.Sum = 0.0;
    helicopter441plots_B.Gain_dg = 0.0;
    helicopter441plots_B.Sum_k = 0.0;
    helicopter441plots_B.Sum2 = 0.0;
    helicopter441plots_B.K_ei = 0.0;
    helicopter441plots_B.Gain_l = 0.0;
    helicopter441plots_B.BackmotorSaturation = 0.0;
    helicopter441plots_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter441plots_M->ModelData.defaultParam = ((real_T *)
    &helicopter441plots_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter441plots_X;
    helicopter441plots_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter441plots_X, 0,
                  sizeof(X_helicopter441plots_T));
  }

  /* states (dwork) */
  helicopter441plots_M->ModelData.dwork = ((void *) &helicopter441plots_DW);
  (void) memset((void *)&helicopter441plots_DW, 0,
                sizeof(DW_helicopter441plots_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter441plots_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter441plots_DW.TimeStampA = 0.0;
  helicopter441plots_DW.LastUAtTimeA = 0.0;
  helicopter441plots_DW.TimeStampB = 0.0;
  helicopter441plots_DW.LastUAtTimeB = 0.0;
  helicopter441plots_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter441plots_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter441plots_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter441plots_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter441plots_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter441plots_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter441plots_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter441plots_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter441plots_M->Sizes.numBlocks = (58);/* Number of blocks */
  helicopter441plots_M->Sizes.numBlockIO = (19);/* Number of block outputs */
  helicopter441plots_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helicopter441plots_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
