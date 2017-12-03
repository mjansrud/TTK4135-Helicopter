/*
 * helicopter44.c
 *
 * Code generation for model "helicopter44".
 *
 * Model version              : 1.173
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Apr 05 08:52:16 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter44.h"
#include "helicopter44_private.h"
#include "helicopter44_dt.h"

/* Block signals (auto storage) */
B_helicopter44_T helicopter44_B;

/* Continuous states */
X_helicopter44_T helicopter44_X;

/* Block states (auto storage) */
DW_helicopter44_T helicopter44_DW;

/* Real-time model */
RT_MODEL_helicopter44_T helicopter44_M_;
RT_MODEL_helicopter44_T *const helicopter44_M = &helicopter44_M_;

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
  helicopter44_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter44_output(void)
{
  /* local block i/o variables */
  real_T rtb_Derivative;
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* set solver stop time */
    if (!(helicopter44_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter44_M->solverInfo,
                            ((helicopter44_M->Timing.clockTickH0 + 1) *
        helicopter44_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter44_M->solverInfo,
                            ((helicopter44_M->Timing.clockTick0 + 1) *
        helicopter44_M->Timing.stepSize0 + helicopter44_M->Timing.clockTickH0 *
        helicopter44_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter44_M)) {
    helicopter44_M->Timing.t[0] = rtsiGetT(&helicopter44_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter44/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter44_DW.HILReadEncoderTimebase_Task,
        1, &helicopter44_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter44_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helicopter44_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter44_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter44_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter44_M->Timing.t[0];

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

    helicopter44_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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

  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter44_B.TravelCounttorad = helicopter44_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter44_B.Gain = helicopter44_P.Gain_Gain *
      helicopter44_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    helicopter44_B.Sum1 = helicopter44_P.travel_offsetdeg_Value +
      helicopter44_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.TravelTransferFcn_C *
    helicopter44_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.TravelTransferFcn_D *
    helicopter44_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter44_B.Gain_d = helicopter44_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter44_B.PitchCounttorad = helicopter44_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter44_B.Gain_i = helicopter44_P.Gain_Gain_a *
      helicopter44_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.PitchTransferFcn_C *
    helicopter44_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.PitchTransferFcn_D *
    helicopter44_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter44_B.Gain_b = helicopter44_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter44_B.ElevationCounttorad = helicopter44_P.ElevationCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter44_B.Gain_e = helicopter44_P.Gain_Gain_lv *
      helicopter44_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter44_B.Sum = helicopter44_B.Gain_e +
      helicopter44_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter44_P.ElevationTransferFcn_C *
    helicopter44_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter44_P.ElevationTransferFcn_D *
    helicopter44_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter44_B.Gain_dg = helicopter44_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helicopter44_P.Gain1_Gain * helicopter44_B.Sum;
  rtb_Gain1_idx_5 = helicopter44_P.Gain1_Gain * helicopter44_B.Gain_dg;

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter44_B.Sum_k = ((rtb_Derivative - helicopter44_P.Gain1_Gain *
    helicopter44_B.Gain_i) * helicopter44_P.K_pp - helicopter44_P.Gain1_Gain *
    helicopter44_B.Gain_b * helicopter44_P.K_pd) + helicopter44_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter44_X.Integrator_CSTATE >= helicopter44_P.Integrator_UpperSat ) {
    helicopter44_X.Integrator_CSTATE = helicopter44_P.Integrator_UpperSat;
  } else if (helicopter44_X.Integrator_CSTATE <=
             (helicopter44_P.Integrator_LowerSat) ) {
    helicopter44_X.Integrator_CSTATE = (helicopter44_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter44_X.Integrator_CSTATE;

  /* FromWorkspace: '<Root>/From Workspace2' */
  {
    real_T *pDataValues = (real_T *)
      helicopter44_DW.FromWorkspace2_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter44_DW.FromWorkspace2_PWORK.TimePtr;
    int_T currTimeIndex = helicopter44_DW.FromWorkspace2_IWORK.PrevIndex;
    real_T t = helicopter44_M->Timing.t[0];

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

    helicopter44_DW.FromWorkspace2_IWORK.PrevIndex = currTimeIndex;

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
  rtb_Derivative -= rtb_Gain1_idx_4;

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter44_B.Sum2 = ((helicopter44_P.K_ep * rtb_Derivative + rtb_Backgain) -
    helicopter44_P.K_ed * rtb_Gain1_idx_5) + helicopter44_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter44_B.Sum2 - helicopter44_B.Sum_k) *
    helicopter44_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter44_B.K_ei = helicopter44_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter44_DW.TimeStampA >= helicopter44_M->Timing.t[0]) &&
      (helicopter44_DW.TimeStampB >= helicopter44_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Gain1_idx_4 = helicopter44_DW.TimeStampA;
    lastU = &helicopter44_DW.LastUAtTimeA;
    if (helicopter44_DW.TimeStampA < helicopter44_DW.TimeStampB) {
      if (helicopter44_DW.TimeStampB < helicopter44_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter44_DW.TimeStampB;
        lastU = &helicopter44_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter44_DW.TimeStampA >= helicopter44_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter44_DW.TimeStampB;
        lastU = &helicopter44_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter44_B.PitchCounttorad - *lastU) /
      (helicopter44_M->Timing.t[0] - rtb_Gain1_idx_4);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter44_B.Gain_l = helicopter44_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter44_P.BackmotorSaturation_UpperSat) {
    helicopter44_B.BackmotorSaturation =
      helicopter44_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter44_P.BackmotorSaturation_LowerSat) {
    helicopter44_B.BackmotorSaturation =
      helicopter44_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter44_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_4 = (helicopter44_B.Sum_k + helicopter44_B.Sum2) *
    helicopter44_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > helicopter44_P.FrontmotorSaturation_UpperSat) {
    helicopter44_B.FrontmotorSaturation =
      helicopter44_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < helicopter44_P.FrontmotorSaturation_LowerSat) {
    helicopter44_B.FrontmotorSaturation =
      helicopter44_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter44_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter44/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter44_DW.HILWriteAnalog_Buffer[0] =
        helicopter44_B.FrontmotorSaturation;
      helicopter44_DW.HILWriteAnalog_Buffer[1] =
        helicopter44_B.BackmotorSaturation;
      result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILWriteAnalog_channels, 2,
        &helicopter44_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter44_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter44_DW.TimeStampA == (rtInf)) {
    helicopter44_DW.TimeStampA = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeA;
  } else if (helicopter44_DW.TimeStampB == (rtInf)) {
    helicopter44_DW.TimeStampB = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeB;
  } else if (helicopter44_DW.TimeStampA < helicopter44_DW.TimeStampB) {
    helicopter44_DW.TimeStampA = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeA;
  } else {
    helicopter44_DW.TimeStampB = helicopter44_M->Timing.t[0];
    lastU = &helicopter44_DW.LastUAtTimeB;
  }

  *lastU = helicopter44_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter44_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter44_M->solverInfo);
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
  if (!(++helicopter44_M->Timing.clockTick0)) {
    ++helicopter44_M->Timing.clockTickH0;
  }

  helicopter44_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter44_M->solverInfo);

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
    if (!(++helicopter44_M->Timing.clockTick1)) {
      ++helicopter44_M->Timing.clockTickH1;
    }

    helicopter44_M->Timing.t[1] = helicopter44_M->Timing.clockTick1 *
      helicopter44_M->Timing.stepSize1 + helicopter44_M->Timing.clockTickH1 *
      helicopter44_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter44_derivatives(void)
{
  XDot_helicopter44_T *_rtXdot;
  _rtXdot = ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter44_P.TravelTransferFcn_A *
    helicopter44_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter44_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter44_P.PitchTransferFcn_A *
    helicopter44_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter44_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter44_P.ElevationTransferFcn_A *
    helicopter44_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter44_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter44_X.Integrator_CSTATE <=
            (helicopter44_P.Integrator_LowerSat) );
    usat = ( helicopter44_X.Integrator_CSTATE >=
            helicopter44_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter44_B.K_ei > 0)) ||
        (usat && (helicopter44_B.K_ei < 0)) ) {
      ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter44_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter44_T *) helicopter44_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter44_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter44/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter44_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter44_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter44_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
      return;
    }

    if ((helicopter44_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter44_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter44_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter44_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter44_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_input_chan, 8U,
        &helicopter44_DW.HILInitialize_AIMinimums[0],
        &helicopter44_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter44_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter44_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter44_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter44_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_output_cha, 8U,
        &helicopter44_DW.HILInitialize_AOMinimums[0],
        &helicopter44_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_analog_output_cha, 8U,
        &helicopter44_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if (helicopter44_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_analog_output_cha, 8U,
         &helicopter44_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter44_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter44_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter44_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter44_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter44_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter44_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter44_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_encoder_channels, 8U,
        &helicopter44_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter44_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter44_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter44_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter44_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter44_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter44_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter44_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter44_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helicopter44_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter44_DW.HILInitialize_Card,
          &helicopter44_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter44_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter44_DW.HILInitialize_Card,
          &helicopter44_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter44_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter44_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter44_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter44_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter44_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter44_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter44_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter44_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter44_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter44_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter44_DW.HILInitialize_POSortedFreqs
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter44_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        &helicopter44_DW.HILInitialize_POSortedFreqs[0],
        &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter44_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter44_DW.HILInitialize_Card,
        helicopter44_P.HILInitialize_pwm_channels, 8U,
        &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }

    if (helicopter44_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter44_DW.HILInitialize_Card,
         helicopter44_P.HILInitialize_pwm_channels, 8U,
         &helicopter44_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter44/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter44_DW.HILInitialize_Card,
      helicopter44_P.HILReadEncoderTimebase_samples_,
      helicopter44_P.HILReadEncoderTimebase_channels, 3,
      &helicopter44_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter44_M, _rt_error_message);
    }
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

    helicopter44_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter44_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter44_DW.FromWorkspace_IWORK.PrevIndex = 0;
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

    helicopter44_DW.FromWorkspace2_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter44_DW.FromWorkspace2_PWORK.DataPtr = (void *) pDataValues0;
    helicopter44_DW.FromWorkspace2_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter44_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter44_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter44_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter44_X.Integrator_CSTATE = helicopter44_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter44_DW.TimeStampA = (rtInf);
  helicopter44_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter44_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter44/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter44_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter44_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter44_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter44_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter44_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter44_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter44_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter44_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter44_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter44_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter44_DW.HILInitialize_Card
                         , helicopter44_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter44_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter44_DW.HILInitialize_AOVoltages[0]
                         , &helicopter44_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter44_DW.HILInitialize_Card,
            helicopter44_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter44_DW.HILInitialize_AOVoltages
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter44_DW.HILInitialize_Card,
            helicopter44_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter44_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter44_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter44_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter44_DW.HILInitialize_Card);
    hil_close(helicopter44_DW.HILInitialize_Card);
    helicopter44_DW.HILInitialize_Card = NULL;
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
  helicopter44_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter44_update();
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
  helicopter44_initialize();
}

void MdlTerminate(void)
{
  helicopter44_terminate();
}

/* Registration function */
RT_MODEL_helicopter44_T *helicopter44(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter44_P.Integrator_UpperSat = rtInf;
  helicopter44_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter44_M, 0,
                sizeof(RT_MODEL_helicopter44_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter44_M->solverInfo,
                          &helicopter44_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter44_M->solverInfo, &rtmGetTPtr(helicopter44_M));
    rtsiSetStepSizePtr(&helicopter44_M->solverInfo,
                       &helicopter44_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter44_M->solverInfo, &helicopter44_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter44_M->solverInfo, (real_T **)
                         &helicopter44_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter44_M->solverInfo,
      &helicopter44_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter44_M->solverInfo, (&rtmGetErrorStatus
      (helicopter44_M)));
    rtsiSetRTModelPtr(&helicopter44_M->solverInfo, helicopter44_M);
  }

  rtsiSetSimTimeStep(&helicopter44_M->solverInfo, MAJOR_TIME_STEP);
  helicopter44_M->ModelData.intgData.f[0] = helicopter44_M->ModelData.odeF[0];
  helicopter44_M->ModelData.contStates = ((real_T *) &helicopter44_X);
  rtsiSetSolverData(&helicopter44_M->solverInfo, (void *)
                    &helicopter44_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter44_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter44_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter44_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter44_M->Timing.sampleTimes =
      (&helicopter44_M->Timing.sampleTimesArray[0]);
    helicopter44_M->Timing.offsetTimes =
      (&helicopter44_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter44_M->Timing.sampleTimes[0] = (0.0);
    helicopter44_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter44_M->Timing.offsetTimes[0] = (0.0);
    helicopter44_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter44_M, &helicopter44_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter44_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter44_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter44_M, -1);
  helicopter44_M->Timing.stepSize0 = 0.002;
  helicopter44_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter44_M->Sizes.checksums[0] = (377538251U);
  helicopter44_M->Sizes.checksums[1] = (815097936U);
  helicopter44_M->Sizes.checksums[2] = (2225222808U);
  helicopter44_M->Sizes.checksums[3] = (2344118900U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter44_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter44_M->extModeInfo,
      &helicopter44_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter44_M->extModeInfo,
                        helicopter44_M->Sizes.checksums);
    rteiSetTPtr(helicopter44_M->extModeInfo, rtmGetTPtr(helicopter44_M));
  }

  helicopter44_M->solverInfoPtr = (&helicopter44_M->solverInfo);
  helicopter44_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter44_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter44_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter44_M->ModelData.blockIO = ((void *) &helicopter44_B);

  {
    helicopter44_B.TravelCounttorad = 0.0;
    helicopter44_B.Gain = 0.0;
    helicopter44_B.Sum1 = 0.0;
    helicopter44_B.Gain_d = 0.0;
    helicopter44_B.PitchCounttorad = 0.0;
    helicopter44_B.Gain_i = 0.0;
    helicopter44_B.Gain_b = 0.0;
    helicopter44_B.ElevationCounttorad = 0.0;
    helicopter44_B.Gain_e = 0.0;
    helicopter44_B.Sum = 0.0;
    helicopter44_B.Gain_dg = 0.0;
    helicopter44_B.Sum_k = 0.0;
    helicopter44_B.Sum2 = 0.0;
    helicopter44_B.K_ei = 0.0;
    helicopter44_B.Gain_l = 0.0;
    helicopter44_B.BackmotorSaturation = 0.0;
    helicopter44_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter44_M->ModelData.defaultParam = ((real_T *)&helicopter44_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter44_X;
    helicopter44_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter44_X, 0,
                  sizeof(X_helicopter44_T));
  }

  /* states (dwork) */
  helicopter44_M->ModelData.dwork = ((void *) &helicopter44_DW);
  (void) memset((void *)&helicopter44_DW, 0,
                sizeof(DW_helicopter44_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter44_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter44_DW.TimeStampA = 0.0;
  helicopter44_DW.LastUAtTimeA = 0.0;
  helicopter44_DW.TimeStampB = 0.0;
  helicopter44_DW.LastUAtTimeB = 0.0;
  helicopter44_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter44_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter44_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter44_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter44_M->Sizes.numY = (0);    /* Number of model outputs */
  helicopter44_M->Sizes.numU = (0);    /* Number of model inputs */
  helicopter44_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter44_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter44_M->Sizes.numBlocks = (55);/* Number of blocks */
  helicopter44_M->Sizes.numBlockIO = (17);/* Number of block outputs */
  helicopter44_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helicopter44_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
