/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "angle_harmonic_propel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleHarmonicPropel::AngleHarmonicPropel(LAMMPS *_lmp) : Angle(_lmp), k(nullptr), theta0(nullptr)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

AngleHarmonicPropel::~AngleHarmonicPropel()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(fpropel);
  }
}

/* ---------------------------------------------------------------------- */

void AngleHarmonicPropel::compute(int eflag, int vflag)
{
  int i1, i2, i3, n, type;
  double delx1, dely1, delz1, delx2, dely2, delz2, delx, dely, delz;
  double eangle, f1[3], f3[3];
  double dtheta, tk;
  double rsq1, rsq2, rsq, r1, r2, r, c, s, a, a11, a12, a22;

  eangle = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int reversal_flag, reversal_type;
  int reversal_index = atom->find_custom("reversal", reversal_flag, reversal_type);
  int* reversal = (reversal_index != -1) ? atom->ivector[reversal_index] : nullptr;

  int theta_flag, theta_type;
  int theta_index = atom->find_custom("theta", theta_flag, theta_type);
  auto* theta = (theta_index != -1) ? atom->dvector[theta_index] : nullptr;

  int force_flag, force_type;
  int force_index = atom->find_custom(force_array.c_str(), force_flag, force_type);
  auto* f_arr = (force_peratom_flag) ? atom->dvector[force_index] : nullptr;

  int nbond_flag, nbond_type;
  int nbond_index = atom->find_custom(nbond_array.c_str(), nbond_flag, nbond_type);
  auto* nbond_arr = (end_correction) ? atom->ivector[nbond_index] : nullptr;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    r2 = sqrt(rsq2);

    // vector between first and last atom

    delx = x[i1][0] - x[i3][0];
    dely = x[i1][1] - x[i3][1];
    delz = x[i1][2] - x[i3][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    // angle (cos and sin)

    c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    // force & energy

    dtheta = acos(c) - theta0[type];
    tk = k[type] * dtheta;

    if (eflag) eangle = tk * dtheta;

    a = -2.0 * tk * s;
    a11 = a * c / rsq1;
    a12 = -a / (r1 * r2);
    a22 = a * c / rsq2;

    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    // apply force on the middle atom
    // reversal along tangent - default to 0 if no reversal vector

    const int rev = (reversal) ? reversal[i2] : 0;
    const int dir = 1 - 2 * rev; // 1 if rev1=0, -1 if rev1=1
    
    double fp = (force_peratom_flag) ? f_arr[i2] / r : fpropel[type] / r;

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += dir * fp * delx;
      f[i2][1] += dir * fp * dely;
      f[i2][2] += dir * fp * delz;
    }

    // update tan theta
    if (theta) {
      double theta_new = std::atan2(dir*dely, dir*delx);
      theta[i2] = theta_new;
      // if end_correction is not set update all atoms in
      // angle so as not to exclude the head/tail atoms.
      if (!end_correction) {
        theta[i1] = theta_new;
        theta[i3] = theta_new;
      }
    }

    // end correction: There are N atoms per filament and only N-2 angles. If the first or last atom
    // in the angle list has only one bond, it is a head/tail atom; we handle them as special cases.
    if (end_correction) {
      // head bond
      if (nbond_arr[i1] == 1) {
        double fp_head = (force_peratom_flag) ? f_arr[i1] / r1 : fpropel[type] / r1;
        if (newton_bond || i1 < nlocal) {
          f[i1][0] += dir * fp_head * delx1;
          f[i1][1] += dir * fp_head * dely1;
          f[i1][2] += dir * fp_head * delz1;
        }
        if (theta) {
          double theta_new = std::atan2(dir*dely1, dir*delx1);
          theta[i1] = theta_new;
        }
      }
      // tail atom
      if (nbond_arr[i3] == 1) {
        double fp_tail = (force_peratom_flag) ? f_arr[i3] / r2 : fpropel[type] / r2;
        if (newton_bond || i3 < nlocal) {
          f[i3][0] -= dir * fp_tail * delx2;
          f[i3][1] -= dir * fp_tail * dely2;
          f[i3][2] -= dir * fp_tail * delz2;
        }
        if (theta) {
          double theta_new = std::atan2(-dir*dely2, -dir*delx2);
          theta[i3] = theta_new;
        }
      }
    }

    if (evflag)
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, f1, f3, delx1, dely1, delz1, delx2, dely2,
               delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleHarmonicPropel::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(k, np1, "angle:k");
  memory->create(theta0, np1, "angle:theta0");
  memory->create(fpropel, np1, "angle:fpropel");

  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleHarmonicPropel::coeff(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR, "Incorrect args for angle coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double theta0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double fp_one = 0;

  if (utils::strmatch(arg[3], "^d_"))
  { // third argument is an array
    force_peratom_flag = true;
    force_array =  std::string(&arg[3][2]);
  }
  else { // thrid argument is a number
    fp_one = utils::numeric(FLERR, arg[3], false, lmp);
  }

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "end_correction") == 0) {
      end_correction = true;
      if (utils::strmatch(arg[iarg+1], "^i_")) nbond_array =  std::string(&arg[iarg+1][2]);
      else error->all(FLERR, "Expected nbond array for end correction" + utils::errorurl(21));
    }
    iarg++;
  }

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = DEG2RAD * theta0_one;
    fpropel[i] = fp_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients" + utils::errorurl(21));
}

/* ---------------------------------------------------------------------- */

double AngleHarmonicPropel::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
   NOTE: cannot write peratom force array
------------------------------------------------------------------------- */

void AngleHarmonicPropel::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&theta0[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&fpropel[1], sizeof(double), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
   NOTE: cannot read peratom force array
------------------------------------------------------------------------- */

void AngleHarmonicPropel::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &theta0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &fpropel[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta0[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&fpropel[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
   NOTE: cannot write peratom force array
------------------------------------------------------------------------- */

void AngleHarmonicPropel::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, k[i], RAD2DEG * theta0[i], fpropel[i]);
}

/* ---------------------------------------------------------------------- */

double AngleHarmonicPropel::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(FLERR, delx1, dely1, delz1);
  double r1 = sqrt(delx1 * delx1 + dely1 * dely1 + delz1 * delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(FLERR, delx2, dely2, delz2);
  double r2 = sqrt(delx2 * delx2 + dely2 * dely2 + delz2 * delz2);

  double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
  c /= r1 * r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  double tk = k[type] * dtheta;
  return tk * dtheta;
}

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *AngleHarmonicPropel::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "k") == 0) return (void *) k;
  if (strcmp(str, "theta0") == 0) return (void *) theta0;
  if (strcmp(str, "fpropel") == 0) return (void *) fpropel;
  return nullptr;
}
