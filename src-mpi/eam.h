/// \file
/// Compute forces for the Embedded Atom Model (EAM).

#ifndef __EAM_H
#define __EAM_H

#include "mytype.h"

struct BasePotentialSt;
struct LinkCellSt;

/// Pointers to the data that is needed in the load and unload functions
/// for the force halo exchange.
/// \see loadForceBuffer
/// \see unloadForceBuffer
typedef struct ForceExchangeDataSt
{
   real_t* dfEmbed; //<! derivative of embedding energy
   struct LinkCellSt* boxes;
}ForceExchangeData;

// typedef struct InterpolationObjectSt
// {
//    int n;          //!< the number of values in the table
//    real_t x0;      //!< the starting ordinate range
//    real_t invDx;   //!< the inverse of the table spacing
//    real_t* values; //!< the abscissa values
// } InterpolationObject;


struct BasePotentialSt* initEamPot(const char* dir, const char* file, const char* type);
#endif
