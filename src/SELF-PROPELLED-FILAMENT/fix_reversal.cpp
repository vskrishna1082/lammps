#include "fix_reversal.h"

#include "atom.h"
#include "error.h"
#include "update.h"
#include <Random123/uniform.hpp>

using namespace LAMMPS_NS;

FixReversal::FixReversal(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix reversal", error); 

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) {
    error->all(FLERR, "Illegal Fix reversal ", nevery);
  }

  if (narg >= 6) {
    avg_runtime = utils::numeric(FLERR, arg[4], false, lmp);
    uk[0] = utils::numeric(FLERR, arg[5], false, lmp);
    if (narg >= 7)
      tau_peratom_flag = utils::inumeric(FLERR, arg[6], false, lmp); // overrides default tau
  } else  {
    error->all(FLERR, "Illegal Fix reversal command");
  }

  if (avg_runtime < 0.0) {
    error->all(FLERR, "Illegal average runtime for Fix reversal");
  }
  rate = update->dt / avg_runtime;

}

FixReversal::~FixReversal()
{
  // destructor
}

int FixReversal::setmask()
{
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

void FixReversal::init()
{
  if (atom->molecule_flag == 0) error->all(FLERR, "Fix reversal requires molecule IDs");
}

void FixReversal::end_of_step()
{
  int* mask = atom->mask;
  tagint* mol = atom->molecule;
  const int nlocal = atom->nlocal;


  int reversal_flag;
  int reversal_type;
  int reversal_index = atom->find_custom("reversal", reversal_flag, reversal_type);
  int taus_flag;
  int taus_type;
  int tau_index = atom->find_custom("tau", taus_flag, taus_type);
  // reversal vector needs to exist for fix to work
  if (reversal_index == -1) error->all(FLERR, "Fix reversal needs custom reversal variable");

  int* reversal = atom->ivector[reversal_index];
  double* taus = atom->dvector[tau_index];

  // flip value of reversal for all atoms in group
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit) 
    { 
      double rate_i = tau_peratom_flag ? update->dt / taus[i] : rate;
      RNG::ctr_type c = {{0}};
      c[0] = update->ntimestep;
      RNG::key_type k = uk;
      k[0] += mol[i];
      auto r = rng(c,k);
      if (r123::u01<double>(r.v[0]) <= rate_i) {
        reversal[i] = !reversal[i];
      }
    }
  }
}
