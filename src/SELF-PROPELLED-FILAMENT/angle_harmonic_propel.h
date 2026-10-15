/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(harmonic/propel,AngleHarmonicPropel);
// clang-format on
#else

#ifndef LMP_ANGLE_HARMONIC_PROPEL_H
#define LMP_ANGLE_HARMONIC_PROPEL_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleHarmonicPropel : public Angle {
 public:
  AngleHarmonicPropel(class LAMMPS *);
  ~AngleHarmonicPropel() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;
  void *extract(const char *, int &) override;

 protected:
  double *k, *theta0, *fpropel;
  std::string force_array, nbond_array;
  bool force_peratom_flag = false, end_correction = false;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
