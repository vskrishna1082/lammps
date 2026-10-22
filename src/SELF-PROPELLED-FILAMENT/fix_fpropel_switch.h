#ifdef FIX_CLASS
// clang-format off
FixStyle(fpropel/switch,FixFpropelSwitch);
// clang-format on
#else

#ifndef LMP_FIX_FPROPEL_SWITCH_H
#define LMP_FIX_FPROPEL_SWITCH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFpropelSwitch : public Fix {
 public:
  FixFpropelSwitch(class LAMMPS *, int, char **);
  ~FixFpropelSwitch() override;
  int setmask() override;
  void end_of_step() override;
 protected:
  double activation_rate, deactivation_rate, active_force, inactive_force = 0;
  class RanMars *rng;
  int seed = 6975742;
  std::string force_array;
};
}

#endif
#endif
