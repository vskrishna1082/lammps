#include "fix_reversal.h"

#include "atom.h"
#include "error.h"
#include "update.h"
#include <Random123/uniform.hpp>

using namespace LAMMPS_NS;

FixReversal::FixReversal(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix reversal", error); 
  avg_runtime = utils::numeric(FLERR, arg[4], false, lmp); // average runtime for molecule
  uk[0] = utils::numeric(FLERR, arg[5], false, lmp);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) {
    error->all(FLERR, "Illegal Fix reversal ", nevery);
  }

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tau_peratom") == 0) {
      tau_peratom_flag = 1;
      if (utils::strmatch(arg[iarg+1], "^d_")) {
        int tmp1, tmp2;
        tau_index = atom->find_custom(&arg[iarg+1][2], tmp1, tmp2);
      } else {
        error->all(FLERR, "Illegal Fix reversal command: tau_peratom requires a per_atom vector");
      }
    }
    if (strcmp(arg[iarg],"alignment_factor") == 0) {
      alignment_factor_flag = 1;
      if (utils::strmatch(arg[iarg+1], "^d_")) {
        int tmp1, tmp2;
        af_index = atom->find_custom(&arg[iarg+1][2], tmp1, tmp2);
      } else {
        error->all(FLERR, "Illegal Fix reversal command: alignment_factor requires a per_atom vector");
      }
    }
    if (strcmp(arg[iarg],"molecular") == 0) {
      if (tau_peratom_flag) {
        error->warning(FLERR, "Fix reversal: molecular overrides tau_peratom option");
        tau_peratom_flag = 0;
      }
      if (alignment_factor_flag) {
        error->warning(FLERR, "Fix reversal: molecular overrides alignment_factor option");
        alignment_factor_flag = 0;
      }
      reverse_molecular_flag = 1;
    }
    iarg += 1;
  }

  if (avg_runtime < 0.0) {
    error->all(FLERR, "Illegal average runtime for Fix reversal");
  }
  rate = (update->dt * nevery) / avg_runtime;

  int reversal_flag;
  int reversal_type;
  reversal_index = atom->find_custom("reversal", reversal_flag, reversal_type);
  if (reversal_index == -1) error->all(FLERR, "Fix reversal needs custom reversal variable");
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
  // Assumes number of molecules does not change
  int* mol = atom->molecule;
  int nlocal = atom->nlocal;
  // TO-DO: This whole business isn't thread safe, and has to be reritten!
  molid_idx.clear();
  int nlocalmol = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mol[i] > 0) {
      if (molid_idx.find(mol[i]) == molid_idx.end()) {
        molid_idx[mol[i]] = nlocalmol;
        nlocalmol++;
      }
    }
  }
  std::cout << "NLOCALMOL is " << nlocalmol << std::endl;
  MPI_Allreduce(&nlocalmol, &nmol, 1, MPI_INT, MPI_SUM, world);
  std::cout << "Fix Reversal: Found " << nmol << " Molecules.\n";
  std::cout << "Fix Reversal: Found " << molid_idx.size() << " Molecules.\n";

  mol_reversed.resize(nmol, false);
}

void FixReversal::end_of_step()
{
  if (reverse_molecular_flag) {
    reverse_molecular();
  } else {
    reverse_atomic();
  }
}

/* Reverse each atom independently, then use an MPI_Allreduce to determine which molecules have at least one atom reversed, and reverse all atoms in those molecules. This is useful when atoms have different tau values, and thus different probabilities of reversal. */

void FixReversal::reverse_atomic()
{
  int* mask = atom->mask;
  int* mol = atom->molecule;
  const int nlocal = atom->nlocal;

  int* reversal = atom->ivector[reversal_index];
  double* taus = (tau_peratom_flag) ? atom->dvector[tau_index] : nullptr;
  double* afs = (alignment_factor_flag) ? atom->dvector[af_index] : nullptr;

  std::vector<int> one_mol_reversed(nmol, 0);

  // flip value of reversal for all atoms in group
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    { 
      double rate_i = (taus) ? nevery * update->dt / taus[i] : rate;
      rate_i = (afs) ? rate_i * (1 -  afs[i]) : rate_i;
      RNG::ctr_type c = {{0}};
      c[0] = update->ntimestep;
      RNG::key_type k = uk;
      k[0] += i;
      auto r = rng(c,k);
      // TODO: Correctly calculate the size of each molecule
      if (r123::u01<double>(r.v[0]) <= rate_i / 50) {
        one_mol_reversed[molid_idx[mol[i]]] = 1;
      }
    }
  }
  for (int i = 0; i < nmol; i++) {
    MPI_Allreduce(&one_mol_reversed[i], &mol_reversed[i], 1, MPI_INT, MPI_LOR, world);
  }
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
      if (mol_reversed[molid_idx[mol[i]]]) {
        reversal[i] = !reversal[i];
        if (alignment_factor_flag) {
          afs[i] = -afs[i];
        }
      }
    }
  }
}

/* Reverse the direction at once for each molecule
*  using a single random number to decide whether to reverse each molecule or not. */

void FixReversal::reverse_molecular()
{
  int* mask = atom->mask;
  int* mol = atom->molecule;
  const int nlocal = atom->nlocal;

  int* reversal = atom->ivector[reversal_index];

  // flip value of reversal for all atoms in group
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    { 
      RNG::ctr_type c = {{0}};
      c[0] = update->ntimestep;
      RNG::key_type k = uk;
      k[0] += mol[i];
      auto r = rng(c,k);
      if (r123::u01<double>(r.v[0]) <= rate) {
        reversal[i] = !reversal[i];
      }
    }
  }
}
