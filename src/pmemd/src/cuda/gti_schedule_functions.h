
#ifndef __GTI_SCHEDULE_H__
#define __GTI_SCHEDULE_H__
#ifdef GTI

#include "gti_const.cuh"

template <typename T> __TARGET T smooth_tanh_func(T  x, T x0, T scale, unsigned direction) {
  T ss = HalfF * Tanh((x - x0) / scale);
  return (direction == 0) ? HalfF + ss : HalfF - ss;
}

template <typename T> __TARGET T d_smooth_tanh_func(T  x, T x0, T scale, unsigned direction) {
  T ss = Tanh((x - x0) / scale);
  T tt = HalfF * (OneF - ss * ss) / scale;
  return (direction == 0) ? tt : -tt;
}

template <typename T> __TARGET T smooth_step1_func(T  x) {
  return (x * x * (ThreeF -TwoF * x ));
}
template <typename T> __TARGET T d_smooth_step1_func(T  x) {
  return x* (SixF - SixF * x);
}

template <typename T> __TARGET T smooth_step2_func(T  x) {
//  return (x * x * x * (TenF - FifteenF * x + SixF * x * x));
  return (x * x * x * Fma(x, (Fma(x, SixF, -FifteenF) ), TenF) );
}
template <typename T> __TARGET T d_smooth_step2_func(T  x) {
  return (x * x * (30 - 60 * x + 30 * x * x));
}

template <typename T> __TARGET T smooth_step3_func(T  x) {
  return (x * x * x * x * (35 - 84 * x + 70 * x * x - 20 * x * x *x));
}
template <typename T> __TARGET T d_smooth_step3_func(T  x) {
  return (x * x * x * (140 - 420 * x + 420 * x * x - 140 * x * x * x));
}

template <typename T> __TARGET T smooth_step4_func(T  x) {
  return (x * x * x * x * x* (126 - 420 * x + 540 * x * x - 315 * x * x * x + 70 *x*x*x*x ));
}
template <typename T> __TARGET T d_smooth_step4_func(T  x) {
  return (x * x * x * x* (126*5 - 420 *6* x + 540 *7* x * x - 315 * 8* x * x * x + 70 * 9* x * x * x * x));
}

template <typename T> __TARGET T smooth_step_func(T  x, int n) {
  switch (n) {
  case(-1): return smooth_tanh_func(x, T(0.5), T(0.15), unsigned(One));
  case(0): return  x;
  case(1): return (smooth_step1_func(x));
  case(2): return (smooth_step2_func(x));
  case(3): return (smooth_step3_func(x));
  case(4): return (smooth_step4_func(x));
  }
  
  // Return something to turn off warning
  return x;
}

template <typename T> __TARGET T d_smooth_step_func(T  x, int n) {
  switch (n) {
  case(-1): return d_smooth_tanh_func(x, T(0.5), T(0.15), unsigned(One));
  case(0): return OneF;
  case(1): return (d_smooth_step1_func(x));
  case(2): return (d_smooth_step2_func(x));
  case(3): return (d_smooth_step3_func(x));
  case(4): return (d_smooth_step4_func(x));
  }
  
  // Return something to turn off warning
  return OneF;
}
// Singleton pattern for the global instance 

class Schedule {

public:

  // This section much be consistent w/ gti.F90 and ti.F90
  enum InteractionType { TypeGen = 0, TypeBAT = 1, TypeEleRec = 2, TypeEleCC = 3, TypeEleSC = 4, TypeVDW = 5, TypeRestBA = 6, TypeEleSSC=7, TypeRMSD=8, TypeTotal = 9 };
  enum FunctionType { linear = 0, smooth_step1 = 1, smooth_step2 = 2, smooth_step3 = 3, smooth_step4 = 4, smooth_tanh = -1 };
  enum ScheduleType : unsigned {Lambda=0, Tau=1};
  enum Direction { forward=0, backward=1, external=0, internal=1 };
  enum Dependence { normal=0, reversed=1 };

  void virtual Setup(int interaction, FunctionType type, unsigned match, double p0, double p1)=0;

  void virtual Init(char* schFileName);
  void virtual Init(bool useSchedule);

  void GetWeight(InteractionType interaction, double lambda, double& weight, double& dWeight, Direction direction = forward) const;
  void GetWeight(InteractionType interaction, int sizeN, double lambdas[], double weights[], double dWeights[], Direction direction = forward) const;

protected:

  Schedule() {};
  virtual ~Schedule() {};

  ScheduleType GetScheduleType() const;

  struct Interaction {
    FunctionType functionType = linear;
    Dependence dependence = normal;
    double parameter[4] = { 0.0, 1.0, 0.0, 0.0 };
  };

  bool m_init = false;
};



class LambdaSchedule: public Schedule {

public:

  enum Matchtype: unsigned { symmetric = 0, complementary = 1 };
  FunctionType GetFunctionType(InteractionType interaction) const;

  static LambdaSchedule& GetReference()
  {
    static LambdaSchedule thisInstance;
    return thisInstance;
  };

//protected:

  LambdaSchedule() {};
  virtual ~LambdaSchedule() {};

  void virtual Setup(int interaction, FunctionType type, unsigned match, double p0, double p1);

  struct LambdaInteraction : Interaction {
    Matchtype match = complementary;
  };

  LambdaInteraction m_interaction[TypeTotal];
};



class TauSchedule : public Schedule {

public:

  enum Regiontype : unsigned { external = 0, internal = 1 };
  FunctionType GetFunctionType(InteractionType interaction, Regiontype type = external) const;

  static TauSchedule& GetReference()
  {
    static TauSchedule thisInstance;
    return thisInstance;
  };


//protected:

  TauSchedule() {};
  virtual ~TauSchedule() {};

  void virtual Setup(int interaction, FunctionType type, unsigned match, double p0, double p1);

  struct TauInteraction : Interaction {
  };

  TauInteraction m_interaction[TypeTotal][2];

};


#endif  /* GTI */
#endif  /*  __GTI_SCHEDULE_H__ */

