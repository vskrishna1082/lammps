#include "fix_fpropel_switch.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "random_mars.h"
#include "update.h"

using namespace LAMMPS_NS;

FixFpropelSwitch::FixFpropelSwitch(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix fpropel/switch", error); 

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) {
    error->all(FLERR, "Illegal Fix reversal ", nevery);
  }
  activation_rate = utils::numeric(FLERR, arg[4], false, lmp);
  deactivation_rate = utils::numeric(FLERR, arg[5], false, lmp);
  active_force = utils::numeric(FLERR, arg[6], false, lmp);
  force_array =  std::string(&arg[7][2]);

  // keywords
  int iarg = 8;
  while (iarg < narg) {
    if (utils::strsame(arg[iarg], "seed")) {
      seed = utils::numeric(FLERR, arg[iarg+1], false, lmp);
    }
    if (utils::strsame(arg[iarg], "inactive_force")) {
      inactive_force = utils::numeric(FLERR, arg[iarg+1], false, lmp);
    }
    iarg += 1;
  }
  // initialize Marsaglia RNG with processor-unique seed
  rng = new RanMars(lmp, seed + comm->me);

}

FixFpropelSwitch::~FixFpropelSwitch()
{
  delete rng;
}

int FixFpropelSwitch::setmask()
{
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

void FixFpropelSwitch::end_of_step()
{
  int* mask = atom->mask;
  const int nlocal = atom->nlocal;

  int farr_flag, farr_type;
  int farr_index = atom->find_custom(force_array.c_str(), farr_flag, farr_type);
  if (farr_index == -1) error->all(FLERR, "fix fpropel/switch: No array found with name '" + force_array + "'!");

  auto* force_arr = atom->dvector[farr_index];

  double activation_prob = (update->dt * nevery) * activation_rate;
  double deactivation_prob = (update->dt * nevery) * deactivation_rate;

  // flip force values
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    { 
      bool is_active = (force_arr[i] == active_force) ? true : false;
      double prob = (is_active) ? deactivation_prob : activation_prob;
      if (rng->uniform() <= prob) {
        force_arr[i] = (is_active) ? inactive_force : active_force;
      }
    }
  }
}
