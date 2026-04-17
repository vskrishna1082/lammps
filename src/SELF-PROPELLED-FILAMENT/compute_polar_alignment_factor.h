#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(polar_alignment/atom,ComputePolarAlignment);
// clang-format on
#else

#ifndef LMP_COMPUTE_POLAR_ALIGNMENT_H
#define LMP_COMPUTE_POLAR_ALIGNMENT_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePolarAlignment : public Compute {
 public:
  ComputePolarAlignment(class LAMMPS *, int, char **);
  ~ComputePolarAlignment() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  int nmax, jgroupbit;
  double cutoff_user, cutsq, mycutneigh;
  int cutflag;
  double strength;
  class NeighList *list;
  double *cvec;

};

}    // namespace LAMMPS_NS

#endif
#endif
