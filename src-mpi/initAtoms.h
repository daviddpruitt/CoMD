/// \file
/// Initialize the atom configuration.

#ifndef __INIT_ATOMS_H
#define __INIT_ATOMS_H

#include "mytype.h"

struct SimFlatSt;
struct LinkCellSt;
#ifdef NEIGHBOR_LIST
#define MAX_NEIGHBORS 128  
typedef struct neighborListSt
{
   int neighborList[MAX_NEIGHBORS]; //!< list of neighbors
   unsigned char neighborScale[MAX_NEIGHBORS]; //!< scaling factor for local remote computation
//   real_t distance[MAX_NEIGHBORS];  //!< distance of neighbor
//   real3 distances[MAX_NEIGHBORS]; //!< distances in all dimensions
} neighborList;
#endif

/// Atom data
typedef struct AtomsSt
{
   // atom-specific data
   int nLocal;    //!< total number of atoms on this processor
   int nGlobal;   //!< total number of atoms in simulation

   int* gid;      //!< A globally unique id for each atom
   int* iSpecies; //!< the species index of the atom
 
#ifdef NEIGHBOR_LIST  
   neighborList *neighbors; //!< neighbor information
   int *numNeighbors;       //!< number of neighbors
#endif
   real3*  r;     //!< positions
   real3*  p;     //!< momenta of atoms
   real3*  f;     //!< forces 
   real_t* U;     //!< potential energy per atom
} Atoms;


/// Allocates memory to store atom data.
Atoms* initAtoms(struct LinkCellSt* boxes);
void destroyAtoms(struct AtomsSt* atoms);

void createFccLattice(int nx, int ny, int nz, real_t lat, struct SimFlatSt* s);

void setVcm(struct SimFlatSt* s, real_t vcm[3]);
void setTemperature(struct SimFlatSt* s, real_t temperature);
void randomDisplacements(struct SimFlatSt* s, real_t delta);
#endif
