#include "compute_polar_alignment_factor.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

ComputePolarAlignment::ComputePolarAlignment(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), cvec(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal compute polar_alignment/atom command");

  // doesn't do much; but maybe we might want to group-limit which neighbors
  // are considered
  jgroupbit = group->get_bitmask_by_id(FLERR, "all", "compute coord/atom");

  strength = utils::numeric(FLERR, arg[3], false, lmp);
  double cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  cutsq = cutoff * cutoff;

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
};

/* ---------------------------------------------------------------------- */

ComputePolarAlignment::~ComputePolarAlignment()
{
  if (copymode) return;
  memory->destroy(cvec);
}

/* ---------------------------------------------------------------------- */

void ComputePolarAlignment::init()
{
  if (force->pair == nullptr) {
	error->all(FLERR, "Compute polar_alignment/atom requires a pair style be defined");
  }
  if (sqrt(cutsq) > force->pair->cutforce) {
	error->all(FLERR, "Compute coord/atom cutoff is longer than pairwise cutoff");
  }

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputePolarAlignment::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePolarAlignment::compute_peratom()
{

  invoked_peratom = update->ntimestep;

  if (atom->nmax > nmax) {
	memory->destroy(cvec);
	nmax = atom->nmax;
	memory->create(cvec, nmax, "polar_alignment/atom:cvec");
	vector_atom = cvec;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // compute polar alignment factor for each atom in group
  // dot product of atom orientation and local polar order P
  // P = [sum(cos(n.theta),sum(sin(n.theta)) for n in atom.neighbor]
  // factor = strength*cos(atom.theta)*P[0] + sin(atom.theta)*P[1]

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *mol = atom->molecule;

  // get theta vector
  int theta_flag, theta_type;
  int theta_index = atom->find_custom("theta", theta_flag, theta_type);
  auto *theta = (theta_index != -1) ? atom->dvector[theta_index] : nullptr;

  if (theta == nullptr) {
	error->all(FLERR, "Compute polar_alignment/atom requires a per-atom vector named theta");
  }

  for (int ii = 0; ii < inum; ii++) {
	int i = ilist[ii];
	if (mask[i] & groupbit) {
	  double xtmp = x[i][0];
	  double ytmp = x[i][1];
	  double ztmp = x[i][2];
	  int imol = mol[i];
	  int jnum = numneigh[i];
	  int *jlist = firstneigh[i];

	  double P_x = 0;
	  double P_y = 0;
	  for (int jj = 0; jj < jnum; jj++) {
		int j = jlist[jj];
		j &= NEIGHMASK;

		if (mask[j] & jgroupbit) {
		  if (mol[j] != imol) {
			double delx = xtmp - x[j][0];
			double dely = ytmp - x[j][1];
			double delz = ztmp - x[j][2];
			double rsq = delx * delx + dely * dely + delz * delz;
			if (rsq < cutsq) {
			  P_x += cos(theta[j]);
			  P_y += sin(theta[j]);
			}
		  }
		}
	  }
	  double factor = strength * (cos(theta[i]) * P_x + sin(theta[i]) * P_y);
	  cvec[i] = factor;
	} else {
	  cvec[i] = 0.0;
	}
  }
}

/* ----------------------------------------------------------------------
    memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePolarAlignment::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
