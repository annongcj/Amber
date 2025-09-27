#ifdef GTI
#include <algorithm>    // std::max..
#include <math.h>       // fmax..
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>   
#include <typeinfo>
#include "gti_schedule_functions.h"

using namespace std;

Schedule::ScheduleType Schedule::GetScheduleType() const {

  const std::type_info& myInfo = typeid(*this);
  std::string nameString(myInfo.name());

  if (nameString.find("Lambda") != std::string::npos)
    return Schedule::Lambda;
  else 
    return Schedule::Tau;

}

void Schedule::Init(bool useSchedule) {

  if (!m_init) {
    FunctionType fType;

    const std::type_info&  myInfo = typeid(*this);
    std::string nameString(myInfo.name());

    if (GetScheduleType()==Lambda){
      fType = (useSchedule) ? smooth_step2 : linear;
      for (InteractionType type = TypeGen; type < TypeTotal; type = InteractionType(type + 1)) {
        Setup(type, fType, LambdaSchedule::complementary, ZeroF, OneF);
      }
    } else {
      fType = linear;
      for (InteractionType type = TypeGen; type < TypeTotal; type = InteractionType(type + 1)) {
        Setup(type, fType, TauSchedule::external, ZeroF, OneF); 
        Setup(type, fType, TauSchedule::internal, ZeroF, OneF);
      }
    }

    m_init = true;
  }


};

void Schedule::Init(char* schFileName) {

  if (m_init) return;

  ScheduleType myType = GetScheduleType();

  std::ifstream t;
  std::string fileName(schFileName);
  fileName.erase(fileName.find_last_not_of(" \n\r\t") + 1);

  t.open(fileName.c_str());
  std::vector<std::string> scheule;

  while (t) {
    std::string line;
    std::getline(t, line);
    scheule.push_back(line);
  }
  t.close();

  Schedule::Init(true);  // Load the default first

  if (scheule.size() == 0) {
    m_init = true;
    return;
  }

  bool setSC = false, setSSC = false;

  InteractionType typeSC = TypeGen;
  FunctionType funSC = linear;
  unsigned matchSC= (myType==Lambda) ? LambdaSchedule::complementary : TauSchedule::external;
  double pSC[2] = { ZeroF, OneF };

  for (unsigned i = 0; i < scheule.size(); i++) {
    std::istringstream iss(scheule[i]);
    std::vector<std::string> tokens((std::istream_iterator<std::string>{iss}),
      std::istream_iterator<std::string>());

    InteractionType type = TypeGen;
    FunctionType fun = linear;
    unsigned match = (myType == Lambda) ? LambdaSchedule::complementary : TauSchedule::external;
    double p[2] = { ZeroF, OneF };

    if (tokens.size() >= 3) {
      if (tokens[0].find("Gen") != string::npos) {
        type = TypeGen;
      }
      else if (tokens[0].find("BAT") != string::npos) {
        type = TypeBAT;
      }
      else if (tokens[0].find("EleRec") != string::npos) {
        type = TypeEleRec;
      }
      else if (tokens[0].find("EleCC") != string::npos) {
        type = TypeEleCC;
      }
      else if (tokens[0].find("EleSC") != string::npos) {
        type = TypeEleSC;
        setSC = true;
      }
      else if (tokens[0].find("EleSSC") != string::npos) {
        type = TypeEleSSC;
      }
      else if (tokens[0].find("VDW") != string::npos) {
        type = TypeVDW;
      }
      else if (tokens[0].find("RestBA") != string::npos) {
        type = TypeRestBA;
      }
      else if (tokens[0].find("RMSD") != string::npos) {
        type = TypeRMSD;
      }

      if (tokens[1].find("linear") != string::npos) {
        fun = linear;
      }
      else if (tokens[1].find("smooth_step0") != string::npos) {
        fun = linear;
      }
      else if (tokens[1].find("smooth_tanh") != string::npos) {
        fun = smooth_tanh; p[0] = 0.5; p[1] = 0.2;
      }
      else if (tokens[1].find("smooth_step1") != string::npos) {
        fun = smooth_step1;
      }
      else if (tokens[1].find("smooth_step2") != string::npos) {
        fun = smooth_step2;
      }
      else if (tokens[1].find("smooth_step3") != string::npos) {
        fun = smooth_step3;
      }
      else if (tokens[1].find("smooth_step4") != string::npos) {
        fun = smooth_step4;
      }

      if (myType == Lambda) {
        if (tokens[2].find("symmetric") != string::npos) {
          match = LambdaSchedule::symmetric;
        } else if (tokens[2].find("complementary") != string::npos) {
          match = LambdaSchedule::complementary;
        }
      } else {
        if (tokens[2].find("external") != string::npos) {
          match = TauSchedule::external;
        } else if (tokens[2].find("internal") != string::npos) {
          match = TauSchedule::internal;
        }
      }

      if (tokens.size() >= 5) {
        for (unsigned j = 3; j < 5; j++) {
          std::string& s = tokens[j];
          s.erase(remove(s.begin(), s.end(), ' '), s.end());
          s.erase(remove(s.begin(), s.end(), ','), s.end());
          p[j - 3] = atof(tokens[j].c_str());
        }
      }

      if (type == TypeEleSC && setSC) {
        pSC[0] = p[0]; pSC[1] = p[1];
        typeSC = type;
        funSC = fun;
        matchSC = match;
      }
      //for debug
      printf("Lambda Scheduling: %4d %4d %4d %6.4f %6.4f \n", type, fun, match, p[0], p[1]);
      Setup(type, fun, match, p[0], p[1]);
    }
  }

  if (setSC && !setSSC) {
    // if SC is set and SSC is not set, then SSC will be set to the same as SC  
     Setup(TypeEleSSC, funSC, matchSC, pSC[0], pSC[1]);
  }

  m_init = true;
}

void Schedule::GetWeight(Schedule::InteractionType interaction, double lambda, double& weight, double& dWeight, Direction direction) const {

  double lambdas[1] = { lambda };
  double w[1], dw[1];

  GetWeight(interaction, 1, lambdas, w, dw, direction);
  weight = w[0]; dWeight = dw[0];

}

void Schedule::GetWeight(Schedule::InteractionType interaction, int sizeN, double lambdas[], double weights[], double dWeights[], Direction direction) const {

  Interaction* pInt = NULL;
  
  bool compLambda = false;
  bool compWeight = false;

  if (GetScheduleType() == Lambda) {
    pInt = &((LambdaSchedule*)this)->m_interaction[interaction];
    compLambda = (direction == backward && ((LambdaSchedule::LambdaInteraction*)pInt)->match == LambdaSchedule::symmetric);
    bool b1= (direction == backward && ((LambdaSchedule::LambdaInteraction*)pInt)->match == LambdaSchedule::complementary);
    bool b2 = (pInt->dependence == reversed);
    compWeight =(b1 != b2);
  } else {
    pInt = &((TauSchedule*)this)->m_interaction[interaction][(direction==external) ? 0 : 1];
    compLambda = false;
    compWeight = false;
  }

  if (!pInt) return;

  double ratio = OneF / (pInt->parameter[1] - pInt->parameter[0]);
  double tt;

  std::vector<double> ll(sizeN, Zero);
  std::vector<double> ww(sizeN, Zero);  
  std::vector<double> dw(sizeN, Zero);
   
  for (int i = 0; i < sizeN; i++) {
    ll[i] = (compLambda) ? OneF - lambdas[i]: lambdas[i] ;
  }
  
  switch (pInt->functionType) {
    case(linear):case(smooth_step1):case(smooth_step2):case(smooth_step3):case(smooth_step4):
      for (int i = 0; i < sizeN; i++) {
        if (ll[i] < pInt->parameter[0]) {
          ww[i] = OneF;
          dw[i] = ZeroF;
        } else if (ll[i]> pInt->parameter[1]) {
          ww[i] = ZeroF;
          dw[i] = ZeroF;
        } else {
          tt = (ll[i] - pInt->parameter[0]) * ratio;
          if (pInt->functionType >= linear && pInt->functionType <= smooth_step4) {
            ww[i] = OneF - smooth_step_func(tt, int(pInt->functionType));
            dw[i] = - ratio * d_smooth_step_func(tt, int(pInt->functionType));
          } 
        }
      }
    break;
    case(smooth_tanh):
      for (int i = 0; i < sizeN; i++) {
        ww[i] = smooth_tanh_func<double>(ll[i], pInt->parameter[0], pInt->parameter[1], 1);
        ww[i] = ww[i] * pInt->parameter[2] + pInt->parameter[3];
        dw[i] = d_smooth_tanh_func<double>(ll[i], pInt->parameter[0], pInt->parameter[1], 1);
        dw[i] *= pInt->parameter[2];
      }
    break;
  }

  for (int i = 0; i < sizeN; i++) {
    weights[i] = compWeight ?  OneF - ww[i]: ww[i];
    double tt = (direction == backward) ? -dw[i] : dw[i];
    if (pInt->dependence == reversed) tt = -tt;
    dWeights[i] = tt;
  }
};


void LambdaSchedule::Setup(int interaction, FunctionType type, unsigned match, double p0, double p1) {

  interaction = max(min(interaction, TypeTotal - 1), int(Zero));
  p0 = fmax(fmin(1.0, p0), ZeroF);
  p1 = fmax(fmin(1.0, p1), ZeroF);

  double t0 = p0, t1 = p1, t2, t3;
  if (type == linear || type == smooth_step1 || type == smooth_step2 || type == smooth_step3 || type == smooth_step4) {
    t0 = fmin(p0, p1);
    t1 = fmax(p0, p1);
  }

  m_interaction[interaction].functionType = type;
  m_interaction[interaction].match = (LambdaSchedule::Matchtype)match;
  m_interaction[interaction].dependence = (p0 < p1) ? normal : reversed;
  m_interaction[interaction].parameter[0] = t0;
  m_interaction[interaction].parameter[1] = t1;

  if (type == smooth_tanh) {
    t2 = smooth_tanh_func<double>(ZeroF, t0, t1, 1);
    t3 = smooth_tanh_func<double>(OneF, t0, t1, 1);

    m_interaction[interaction].parameter[2] = OneF / (t2 - t3);
    m_interaction[interaction].parameter[3] = (OneF - t2 / (t2 - t3));

  }

};


Schedule::FunctionType LambdaSchedule::GetFunctionType(Schedule::InteractionType interaction) const {

  return m_interaction[interaction].functionType;
}


void TauSchedule::Setup(int interaction, FunctionType type, unsigned match, double p0, double p1) {

  interaction = max(min(interaction, TypeTotal - 1), int(Zero));
  p0 = fmax(fmin(1.0, p0), ZeroF);
  p1 = fmax(fmin(1.0, p1), ZeroF);

  double t0 = p0, t1 = p1, t2, t3;
  if (type == linear || type == smooth_step1 || type == smooth_step2 || type == smooth_step3 || type == smooth_step4) {
    t0 = fmin(p0, p1);
    t1 = fmax(p0, p1);
  }

  unsigned i = (unsigned) match;
  m_interaction[interaction][i].functionType = type;
  m_interaction[interaction][i].parameter[0] = t0;
  m_interaction[interaction][i].parameter[1] = t1;

  if (type == smooth_tanh) {
    t2 = smooth_tanh_func<double>(ZeroF, t0, t1, 1);
    t3 = smooth_tanh_func<double>(OneF, t0, t1, 1);

    m_interaction[interaction][i].parameter[2] = OneF / (t2 - t3);
    m_interaction[interaction][i].parameter[3] = (OneF - t2 / (t2 - t3));
  }

};

Schedule::FunctionType TauSchedule::GetFunctionType(InteractionType interaction, Regiontype type) const {

  return m_interaction[interaction][(unsigned)type].functionType;

}
#endif /* GTI */
