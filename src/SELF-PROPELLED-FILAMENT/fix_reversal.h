#ifdef FIX_CLASS
// clang-format off
FixStyle(reversal,FixReversal);
// clang-format on
#else

#ifndef LMP_FIX_REVERSAL_H
#define LMP_FIX_REVERSAL_H

#include "fix.h"

#include <Random123/philox.h>
#include <unordered_map>

typedef r123::Philox2x64 RNG;

namespace LAMMPS_NS {

class FixReversal : public Fix {
 public:
  FixReversal(class LAMMPS *, int, char **);
  ~FixReversal() override;
  void init() override;
  int setmask() override;
  void end_of_step() override;
  void reverse_molecular();
  void reverse_atomic();
 protected:
  double avg_runtime;
  double rate;
  int nmol;
  RNG rng;
  RNG::ukey_type uk = {{43923764583}};
  int tau_peratom_flag = 0;
  int alignment_factor_flag = 0;
  int reverse_molecular_flag = 0;
  int reversal_index, tau_index, af_index;
  std::vector<int> mol_reversed;
  std::unordered_map<int,int> molid_idx;
};
}

#endif
#endif
