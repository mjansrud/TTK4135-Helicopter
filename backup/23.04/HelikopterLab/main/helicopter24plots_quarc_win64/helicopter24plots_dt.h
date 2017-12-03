/*
 * helicopter24plots_dt.h
 *
 * Code generation for model "helicopter24plots".
 *
 * Model version              : 1.185
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Apr 05 10:33:53 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "ext_types.h"

/* data type size table */
static uint_T rtDataTypeSizes[] = {
  sizeof(real_T),
  sizeof(real32_T),
  sizeof(int8_T),
  sizeof(uint8_T),
  sizeof(int16_T),
  sizeof(uint16_T),
  sizeof(int32_T),
  sizeof(uint32_T),
  sizeof(boolean_T),
  sizeof(fcn_call_T),
  sizeof(int_T),
  sizeof(pointer_T),
  sizeof(action_T),
  2*sizeof(uint32_T),
  sizeof(t_card),
  sizeof(t_task)
};

/* data type name table */
static const char_T * rtDataTypeNames[] = {
  "real_T",
  "real32_T",
  "int8_T",
  "uint8_T",
  "int16_T",
  "uint16_T",
  "int32_T",
  "uint32_T",
  "boolean_T",
  "fcn_call_T",
  "int_T",
  "pointer_T",
  "action_T",
  "timer_uint32_pair_T",
  "t_card",
  "t_task"
};

/* data type transitions for block I/O structure */
static DataTypeTransition rtBTransitions[] = {
  { (char_T *)(&helicopter24plots_B.FromWorkspace1[0]), 0, 0, 36 }
  ,

  { (char_T *)(&helicopter24plots_DW.HILInitialize_AIMinimums[0]), 0, 0, 70 },

  { (char_T *)(&helicopter24plots_DW.HILInitialize_Card), 14, 0, 1 },

  { (char_T *)(&helicopter24plots_DW.HILReadEncoderTimebase_Task), 15, 0, 1 },

  { (char_T *)(&helicopter24plots_DW.FromWorkspace1_PWORK.TimePtr), 11, 0, 17 },

  { (char_T *)(&helicopter24plots_DW.HILInitialize_ClockModes[0]), 6, 0, 46 },

  { (char_T *)(&helicopter24plots_DW.HILInitialize_POSortedChans[0]), 7, 0, 8 },

  { (char_T *)(&helicopter24plots_DW.FromWorkspace1_IWORK.PrevIndex), 10, 0, 3 }
};

/* data type transition table for block I/O structure */
static DataTypeTransitionTable rtBTransTable = {
  8U,
  rtBTransitions
};

/* data type transitions for Parameters structure */
static DataTypeTransition rtPTransitions[] = {
  { (char_T *)(&helicopter24plots_P.K_ed), 0, 0, 24 },

  { (char_T *)(&helicopter24plots_P.HILReadEncoderTimebase_clock), 6, 0, 10 },

  { (char_T *)(&helicopter24plots_P.HILInitialize_analog_input_chan[0]), 7, 0,
    39 },

  { (char_T *)(&helicopter24plots_P.HILInitialize_active), 8, 0, 35 },

  { (char_T *)(&helicopter24plots_P.TravelCounttorad_Gain), 0, 0, 31 },

  { (char_T *)(&helicopter24plots_P.HILReadEncoderTimebase_Active), 8, 0, 2 }
};

/* data type transition table for Parameters structure */
static DataTypeTransitionTable rtPTransTable = {
  6U,
  rtPTransitions
};

/* [EOF] helicopter24plots_dt.h */
