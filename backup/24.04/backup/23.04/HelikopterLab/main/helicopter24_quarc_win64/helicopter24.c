/*
 * helicopter24.c
 *
 * Code generation for model "helicopter24".
 *
 * Model version              : 1.166
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Mar 15 11:48:02 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter24.h"
#include "helicopter24_private.h"
#include "helicopter24_dt.h"

/* Block signals (auto storage) */
B_helicopter24_T helicopter24_B;

/* Continuous states */
X_helicopter24_T helicopter24_X;

/* Block states (auto storage) */
DW_helicopter24_T helicopter24_DW;

/* Real-time model */
RT_MODEL_helicopter24_T helicopter24_M_;
RT_MODEL_helicopter24_T *const helicopter24_M = &helicopter24_M_;

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
  helicopter24_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter24_output(void)
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
  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* set solver stop time */
    if (!(helicopter24_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter24_M->solverInfo,
                            ((helicopter24_M->Timing.clockTickH0 + 1) *
        helicopter24_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter24_M->solverInfo,
                            ((helicopter24_M->Timing.clockTick0 + 1) *
        helicopter24_M->Timing.stepSize0 + helicopter24_M->Timing.clockTickH0 *
        helicopter24_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter24_M)) {
    helicopter24_M->Timing.t[0] = rtsiGetT(&helicopter24_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter24/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter24_DW.HILReadEncoderTimebase_Task,
        1, &helicopter24_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter24_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter24_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter24_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helicopter24_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter24_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter24_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter24_M->Timing.t[0];

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

    helicopter24_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

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
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter24_B.TravelCounttorad = helicopter24_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter24_B.Gain = helicopter24_P.Gain_Gain *
      helicopter24_B.TravelCounttorad;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24_P.TravelTransferFcn_C *
    helicopter24_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter24_P.TravelTransferFcn_D *
    helicopter24_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter24_B.Gain_d = helicopter24_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter24_B.PitchCounttorad = helicopter24_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter24_B.Gain_i = helicopter24_P.Gain_Gain_a *
      helicopter24_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24_P.PitchTransferFcn_C *
    helicopter24_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter24_P.PitchTransferFcn_D *
    helicopter24_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter24_B.Gain_b = helicopter24_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter24_B.ElevationCounttorad = helicopter24_P.ElevationCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter24_B.Gain_e = helicopter24_P.Gain_Gain_lv *
      helicopter24_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter24_B.Sum = helicopter24_B.Gain_e +
      helicopter24_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter24_P.ElevationTransferFcn_C *
    helicopter24_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter24_P.ElevationTransferFcn_D *
    helicopter24_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter24_B.Gain_dg = helicopter24_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helicopter24_P.Gain1_Gain * helicopter24_B.Sum;
  rtb_Gain1_idx_5 = helicopter24_P.Gain1_Gain * helicopter24_B.Gain_dg;

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter24_B.Sum_k = ((rtb_Derivative - helicopter24_P.Gain1_Gain *
    helicopter24_B.Gain_i) * helicopter24_P.K_pp - helicopter24_P.Gain1_Gain *
    helicopter24_B.Gain_b * helicopter24_P.K_pd) + helicopter24_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter24_X.Integrator_CSTATE >= helicopter24_P.Integrator_UpperSat ) {
    helicopter24_X.Integrator_CSTATE = helicopter24_P.Integrator_UpperSat;
  } else if (helicopter24_X.Integrator_CSTATE <=
             (helicopter24_P.Integrator_LowerSat) ) {
    helicopter24_X.Integrator_CSTATE = (helicopter24_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter24_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter24_P.elevation_ref_Value - rtb_Gain1_idx_4;

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter24_B.Sum2 = ((helicopter24_P.K_ep * rtb_Derivative + rtb_Backgain) -
    helicopter24_P.K_ed * rtb_Gain1_idx_5) + helicopter24_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter24_B.Sum2 - helicopter24_B.Sum_k) *
    helicopter24_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter24_B.K_ei = helicopter24_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter24_DW.TimeStampA >= helicopter24_M->Timing.t[0]) &&
      (helicopter24_DW.TimeStampB >= helicopter24_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Gain1_idx_4 = helicopter24_DW.TimeStampA;
    lastU = &helicopter24_DW.LastUAtTimeA;
    if (helicopter24_DW.TimeStampA < helicopter24_DW.TimeStampB) {
      if (helicopter24_DW.TimeStampB < helicopter24_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter24_DW.TimeStampB;
        lastU = &helicopter24_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter24_DW.TimeStampA >= helicopter24_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helicopter24_DW.TimeStampB;
        lastU = &helicopter24_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter24_B.PitchCounttorad - *lastU) /
      (helicopter24_M->Timing.t[0] - rtb_Gain1_idx_4);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter24_B.Gain_l = helicopter24_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter24_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter24_P.BackmotorSaturation_UpperSat) {
    helicopter24_B.BackmotorSaturation =
      helicopter24_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter24_P.BackmotorSaturation_LowerSat) {
    helicopter24_B.BackmotorSaturation =
      helicopter24_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter24_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter24_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_4 = (helicopter24_B.Sum_k + helicopter24_B.Sum2) *
    helicopter24_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > helicopter24_P.FrontmotorSaturation_UpperSat) {
    helicopter24_B.FrontmotorSaturation =
      helicopter24_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < helicopter24_P.FrontmotorSaturation_LowerSat) {
    helicopter24_B.FrontmotorSaturation =
      helicopter24_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter24_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter24_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter24/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter24_DW.HILWriteAnalog_Buffer[0] =
        helicopter24_B.FrontmotorSaturation;
      helicopter24_DW.HILWriteAnalog_Buffer[1] =
        helicopter24_B.BackmotorSaturation;
      result = hil_write_analog(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILWriteAnalog_channels, 2,
        &helicopter24_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter24_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter24_DW.TimeStampA == (rtInf)) {
    helicopter24_DW.TimeStampA = helicopter24_M->Timing.t[0];
    lastU = &helicopter24_DW.LastUAtTimeA;
  } else if (helicopter24_DW.TimeStampB == (rtInf)) {
    helicopter24_DW.TimeStampB = helicopter24_M->Timing.t[0];
    lastU = &helicopter24_DW.LastUAtTimeB;
  } else if (helicopter24_DW.TimeStampA < helicopter24_DW.TimeStampB) {
    helicopter24_DW.TimeStampA = helicopter24_M->Timing.t[0];
    lastU = &helicopter24_DW.LastUAtTimeA;
  } else {
    helicopter24_DW.TimeStampB = helicopter24_M->Timing.t[0];
    lastU = &helicopter24_DW.LastUAtTimeB;
  }

  *lastU = helicopter24_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter24_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter24_M->solverInfo);
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
  if (!(++helicopter24_M->Timing.clockTick0)) {
    ++helicopter24_M->Timing.clockTickH0;
  }

  helicopter24_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter24_M->solverInfo);

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
    if (!(++helicopter24_M->Timing.clockTick1)) {
      ++helicopter24_M->Timing.clockTickH1;
    }

    helicopter24_M->Timing.t[1] = helicopter24_M->Timing.clockTick1 *
      helicopter24_M->Timing.stepSize1 + helicopter24_M->Timing.clockTickH1 *
      helicopter24_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter24_derivatives(void)
{
  XDot_helicopter24_T *_rtXdot;
  _rtXdot = ((XDot_helicopter24_T *) helicopter24_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter24_P.TravelTransferFcn_A *
    helicopter24_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter24_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter24_P.PitchTransferFcn_A *
    helicopter24_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter24_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter24_P.ElevationTransferFcn_A *
    helicopter24_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter24_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter24_X.Integrator_CSTATE <=
            (helicopter24_P.Integrator_LowerSat) );
    usat = ( helicopter24_X.Integrator_CSTATE >=
            helicopter24_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter24_B.K_ei > 0)) ||
        (usat && (helicopter24_B.K_ei < 0)) ) {
      ((XDot_helicopter24_T *) helicopter24_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter24_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter24_T *) helicopter24_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter24_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter24/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter24_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter24_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter24_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24_M, _rt_error_message);
      return;
    }

    if ((helicopter24_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter24_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter24_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter24_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter24_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter24_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_analog_input_chan, 8U,
        &helicopter24_DW.HILInitialize_AIMinimums[0],
        &helicopter24_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter24_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter24_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter24_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter24_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter24_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_analog_output_cha, 8U,
        &helicopter24_DW.HILInitialize_AOMinimums[0],
        &helicopter24_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter24_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter24_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_analog_output_cha, 8U,
        &helicopter24_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if (helicopter24_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter24_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter24_DW.HILInitialize_Card,
         helicopter24_P.HILInitialize_analog_output_cha, 8U,
         &helicopter24_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter24_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter24_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter24_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter24_DW.HILInitialize_Card,
         helicopter24_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter24_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter24_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter24_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter24_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_encoder_channels, 8U,
        &helicopter24_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter24_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter24_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter24_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter24_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter24_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter24_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter24_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter24_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter24_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter24_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter24_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helicopter24_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter24_DW.HILInitialize_Card,
          &helicopter24_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter24_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter24_DW.HILInitialize_Card,
          &helicopter24_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter24_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter24_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter24_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter24_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter24_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter24_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter24_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter24_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter24_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter24_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter24_DW.HILInitialize_POSortedFreqs
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter24_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter24_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_pwm_channels, 8U,
        &helicopter24_DW.HILInitialize_POSortedFreqs[0],
        &helicopter24_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter24_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter24_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter24_DW.HILInitialize_Card,
        helicopter24_P.HILInitialize_pwm_channels, 8U,
        &helicopter24_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }

    if (helicopter24_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter24_DW.HILInitialize_Card,
         helicopter24_P.HILInitialize_pwm_channels, 8U,
         &helicopter24_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter24/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter24_DW.HILInitialize_Card,
      helicopter24_P.HILReadEncoderTimebase_samples_,
      helicopter24_P.HILReadEncoderTimebase_channels, 3,
      &helicopter24_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter24_M, _rt_error_message);
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

    helicopter24_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter24_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter24_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter24_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter24_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter24_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter24_X.Integrator_CSTATE = helicopter24_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter24_DW.TimeStampA = (rtInf);
  helicopter24_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter24_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter24/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter24_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter24_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter24_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter24_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter24_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter24_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter24_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter24_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter24_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter24_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter24_DW.HILInitialize_Card
                         , helicopter24_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter24_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter24_DW.HILInitialize_AOVoltages[0]
                         , &helicopter24_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter24_DW.HILInitialize_Card,
            helicopter24_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter24_DW.HILInitialize_AOVoltages
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter24_DW.HILInitialize_Card,
            helicopter24_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter24_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter24_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter24_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter24_DW.HILInitialize_Card);
    hil_close(helicopter24_DW.HILInitialize_Card);
    helicopter24_DW.HILInitialize_Card = NULL;
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
  helicopter24_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter24_update();
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
  helicopter24_initialize();
}

void MdlTerminate(void)
{
  helicopter24_terminate();
}

/* Registration function */
RT_MODEL_helicopter24_T *helicopter24(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter24_P.Integrator_UpperSat = rtInf;
  helicopter24_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter24_M, 0,
                sizeof(RT_MODEL_helicopter24_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter24_M->solverInfo,
                          &helicopter24_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter24_M->solverInfo, &rtmGetTPtr(helicopter24_M));
    rtsiSetStepSizePtr(&helicopter24_M->solverInfo,
                       &helicopter24_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter24_M->solverInfo, &helicopter24_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter24_M->solverInfo, (real_T **)
                         &helicopter24_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter24_M->solverInfo,
      &helicopter24_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter24_M->solverInfo, (&rtmGetErrorStatus
      (helicopter24_M)));
    rtsiSetRTModelPtr(&helicopter24_M->solverInfo, helicopter24_M);
  }

  rtsiSetSimTimeStep(&helicopter24_M->solverInfo, MAJOR_TIME_STEP);
  helicopter24_M->ModelData.intgData.f[0] = helicopter24_M->ModelData.odeF[0];
  helicopter24_M->ModelData.contStates = ((real_T *) &helicopter24_X);
  rtsiSetSolverData(&helicopter24_M->solverInfo, (void *)
                    &helicopter24_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter24_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter24_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter24_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter24_M->Timing.sampleTimes =
      (&helicopter24_M->Timing.sampleTimesArray[0]);
    helicopter24_M->Timing.offsetTimes =
      (&helicopter24_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter24_M->Timing.sampleTimes[0] = (0.0);
    helicopter24_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter24_M->Timing.offsetTimes[0] = (0.0);
    helicopter24_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter24_M, &helicopter24_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter24_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter24_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter24_M, -1);
  helicopter24_M->Timing.stepSize0 = 0.002;
  helicopter24_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter24_M->Sizes.checksums[0] = (676312908U);
  helicopter24_M->Sizes.checksums[1] = (2504426561U);
  helicopter24_M->Sizes.checksums[2] = (3573356144U);
  helicopter24_M->Sizes.checksums[3] = (984827586U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter24_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter24_M->extModeInfo,
      &helicopter24_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter24_M->extModeInfo,
                        helicopter24_M->Sizes.checksums);
    rteiSetTPtr(helicopter24_M->extModeInfo, rtmGetTPtr(helicopter24_M));
  }

  helicopter24_M->solverInfoPtr = (&helicopter24_M->solverInfo);
  helicopter24_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter24_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter24_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter24_M->ModelData.blockIO = ((void *) &helicopter24_B);

  {
    helicopter24_B.TravelCounttorad = 0.0;
    helicopter24_B.Gain = 0.0;
    helicopter24_B.Gain_d = 0.0;
    helicopter24_B.PitchCounttorad = 0.0;
    helicopter24_B.Gain_i = 0.0;
    helicopter24_B.Gain_b = 0.0;
    helicopter24_B.ElevationCounttorad = 0.0;
    helicopter24_B.Gain_e = 0.0;
    helicopter24_B.Sum = 0.0;
    helicopter24_B.Gain_dg = 0.0;
    helicopter24_B.Sum_k = 0.0;
    helicopter24_B.Sum2 = 0.0;
    helicopter24_B.K_ei = 0.0;
    helicopter24_B.Gain_l = 0.0;
    helicopter24_B.BackmotorSaturation = 0.0;
    helicopter24_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter24_M->ModelData.defaultParam = ((real_T *)&helicopter24_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter24_X;
    helicopter24_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter24_X, 0,
                  sizeof(X_helicopter24_T));
  }

  /* states (dwork) */
  helicopter24_M->ModelData.dwork = ((void *) &helicopter24_DW);
  (void) memset((void *)&helicopter24_DW, 0,
                sizeof(DW_helicopter24_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter24_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter24_DW.TimeStampA = 0.0;
  helicopter24_DW.LastUAtTimeA = 0.0;
  helicopter24_DW.TimeStampB = 0.0;
  helicopter24_DW.LastUAtTimeB = 0.0;
  helicopter24_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter24_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter24_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter24_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter24_M->Sizes.numY = (0);    /* Number of model outputs */
  helicopter24_M->Sizes.numU = (0);    /* Number of model inputs */
  helicopter24_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter24_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter24_M->Sizes.numBlocks = (53);/* Number of blocks */
  helicopter24_M->Sizes.numBlockIO = (16);/* Number of block outputs */
  helicopter24_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helicopter24_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
