/*
 * helicopter441plots.c
 *
 * Code generation for model "helicopter441plots".
 *
 * Model version              : 1.181
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Apr 26 22:59:20 2017
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
    } else if (t >= pTimeValues[79]) {
      currTimeIndex = 78;
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
              pDataValues += 80;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&helicopter441plots_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 80;
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
            pDataValues += 80;
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
    } else if (t >= pTimeValues[79]) {
      currTimeIndex = 78;
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
        pDataValues += 80;
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
    } else if (t >= pTimeValues[79]) {
      currTimeIndex = 78;
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
        pDataValues += 80;
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
      19.75 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1378421413625261, 3.1262155534579983, 3.1033093000299643,
      3.0666274151911783, 3.0144539223941584, 2.9456562771175667,
      2.8595077632935446, 2.7555515879651526, 2.633505110490284,
      2.4931956060320961, 2.334518576064299, 2.1574113214711188,
      1.9618364840294067, 1.7541234930883733, 1.5423954893370344,
      1.3336732751551714, 1.133222483102116, 0.94471151627325012,
      0.77058276422429517, 0.61241316353019215, 0.47119852158277031,
      0.34755720319912214, 0.24185442447598429, 0.1529582783901568,
      0.079464099431143978, 0.019799342367502162, -0.02773773055502602,
      -0.064887788454194489, -0.093328089670144027, -0.11465070782674336,
      -0.13035492400331725, -0.141805295146489, -0.15018945440076378,
      -0.15650177732827505, -0.1615459762043151, -0.16594410466703374,
      -0.17012574663058364, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, -0.015002048909068425, -0.046506351618112118,
      -0.091625013712135411, -0.14672753935514374, -0.20869397118807928,
      -0.27519058110636674, -0.34459405529608839, -0.41582470131356869,
      -0.48818590989947463, -0.56123801783275185, -0.63470811987118858,
      -0.7084290183727211, -0.78229934976684812, -0.83085196376413351,
      -0.84691201500535551, -0.8348888567274525, -0.80180316821222108,
      -0.75404386731546369, -0.6965150081958198, -0.63267840277641207,
      -0.56485856778968746, -0.49456527353459279, -0.42281111489255141,
      -0.35558458434331003, -0.29397671583605128, -0.23865902825456725,
      -0.19014829169011274, -0.14860023159667385, -0.11376120486379818,
      -0.085290472626397343, -0.062816864706295633, -0.045801484572687012,
      -0.033536637017099108, -0.025249291710045134, -0.020176795504160249,
      -0.017592513850874496, -0.016726567854199657, -0.016698213316866054, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10602875205865553,
      0.22266037932317662, 0.31888147181640647, 0.3894436063114417,
      0.43795507377677845, 0.46997264230390068, 0.49051724877547076,
      0.5034310014147434, 0.5114213858602934, 0.51630439857701838,
      0.51925862127063693, 0.52103115488680807, 0.52208728949977667,
      0.34315133236273704, 0.11350630846667278, -0.084975090785043436,
      -0.23383701023337183, -0.33754449835288891, -0.40659200465123835,
      -0.45117274642326882, -0.47932469171725162, -0.49680615715459348,
      -0.50713098870352036, -0.4751314425498539, -0.43542088513455018,
      -0.39096428871725047, -0.34285535866145234, -0.29364582057204508,
      -0.24622893511512536, -0.20122026182056985, -0.15883487758685549,
      -0.12025820821578632, -0.086683258292976018, -0.058571791500755234,
      -0.035850465879172823, -0.018264715728285297, -0.0061201755796251095,
      -0.00020039904061806054, -0.00037016870632875876, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462212, 0.46652650905808429,
      0.38488436997291953, 0.282248537980141, 0.19404586986134692,
      0.12807027410848892, 0.082178425886280354, 0.051655010557090507,
      0.03196153778219981, 0.01953205086689988, 0.011816890774474417,
      0.0070901344646846073, 0.0042245384518745975, -0.71574382854815843,
      -0.91858009558425713, -0.79392559700686482, -0.59544767779331353,
      -0.41482995247806825, -0.27619002519339769, -0.17832296708812181,
      -0.11260778117593113, -0.069925861749367421, -0.041299326195707656,
      0.12799818461466578, 0.15884222966121489, 0.17782638566919903,
      0.19243572022319252, 0.19683815235762897, 0.18966754182767878,
      0.1800346931782221, 0.16954153693485752, 0.15430667748427665,
      0.13429979969124117, 0.11244586716888313, 0.090885302486329617,
      0.070343000603550115, 0.048578160594640749, 0.023679106156028196,
      -0.00067907866284279287, 0.00052904719012034779, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001556144130800542, 1.4921969912328133E-10,
      0.0071896560518316718, 0.02894435971390593, 0.061510510822098588,
      0.10184533027016013, 0.14749780486205782, 0.19650828925795669,
      0.247323961091358, 0.29872773468420027, 0.34977858972422632,
      0.39976157188725209, 0.44814597978252718, 0.46182554940579834,
      0.44872048395170289, 0.41531319891421886, 0.36688401583691366,
      0.30771039805845324, 0.24123516023970931, 0.17020829516854122,
      0.096806385248972865, 0.036879961015705988, 0.0067850048086874641,
      1.9523545340635239E-20, 0.0016386856649306407, 0.0027331119109765857,
      0.0020540115655257088, 0.0007787850405197818, 2.6780510101477567E-20,
      1.0321408629183986E-5, 0.00084494295036762932, 0.0013628361700970776,
      0.0016043465598761383, 0.0011821178505498173, 0.00037091979499453861,
      4.9949745544664166E-20, 4.2942619513269677E-6, 0.00047475502027454039,
      0.0013880739626466541, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0062245765232021681, -0.0062245759263233717, 0.028758623610447891,
      0.087018814648297027, 0.13026460443277063, 0.16133927779224619,
      0.18260989836759081, 0.19604193758359545, 0.20326268733360531,
      0.20561509437136902, 0.20420342016010423, 0.19993192865210302,
      0.19353763158110043, 0.054718278493084742, -0.052420261816381722,
      -0.13362914014993602, -0.19371673230922079, -0.23669447111384159,
      -0.26590095127497565, -0.28410746028467238, -0.29360763967827341,
      -0.23970569693306754, -0.12037982482807409, -0.027140019234749856,
      0.0065547426597225627, 0.00437770498418378, -0.002716401381803507,
      -0.0051009061000237086, -0.0031151401620791272, 4.1285634516735834E-5,
      0.0033384861669537814, 0.0020715728789177934, 0.00096604155911624273,
      -0.001688914837305283, -0.003244792222221115, -0.0014836791799781542,
      1.7177047805307681E-5, 0.0018818430332928537, 0.0036532757694884553,
      0.0027009756835588493, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      19.75 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      -0.36312474930250771, -0.52359877559829882, -0.52359877559829882,
      -0.52359877559829882, -0.52359877559829882, -0.52359877559829882,
      -0.52359877559829882, -0.52359877559829882, -0.52359877559829882,
      -0.52167869174289472, -0.33368502476664719, -0.32683171561741137,
      -0.27520284662012867, -0.21979979667001606, -0.17171172040205326,
      -0.13299890607884016, -0.094796513170249813, -0.059144901020020089,
      -0.031649244794902082, -0.012082862842641295, 0.001983627833298863,
      0.011638613027321866, 0.017051011406549882, 0.015438250829821904,
      0.0049714448847568236, -0.0098818907381999378, 0.00053658251049876291, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

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
      19.75 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0995932243712347,
      -0.1742881330996, 0.53638903301384733, 1.0471975511965976,
      1.0471975511965976, 1.0471975511965976, 1.0471975511965976,
      1.0471975511965976, 1.0471975511965976, 1.0471975511965976,
      1.0471975511965976, 1.0471975511965976, 1.0471975511965976,
      -1.0471975511965976, -1.0471975511965976, -1.0471975511965976,
      -1.0471975511965976, -1.0471975511965976, -1.0471975511965976,
      -1.0471975511965976, -1.0471975511965976, -0.1417911796212577,
      1.0471975511965976, 1.0471975511965976, 0.43734111818124677,
      -0.008613632169730278, -0.094356196254130845, -0.046284569107760673,
      0.013422642172544174, 0.038821037137737081, 0.052920351057059681,
      -0.006906346532131501, -0.00855726665078601, -0.03725229993618237,
      -0.030045350947998307, 0.016380757637552725, 0.018449902719617314,
      0.029903363959021972, 0.035874590174252367, -0.00014894327664533317, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

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
  helicopter441plots_M->Sizes.checksums[0] = (144244424U);
  helicopter441plots_M->Sizes.checksums[1] = (3442827034U);
  helicopter441plots_M->Sizes.checksums[2] = (2156392845U);
  helicopter441plots_M->Sizes.checksums[3] = (834255105U);

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
