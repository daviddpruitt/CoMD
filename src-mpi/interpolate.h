#include "eam.h"
#include "CoMDTypes.h"
#include <math.h>

/// Interpolate a table to determine f(r) and its derivative f'(r).
///
/// The forces on the particle are much more sensitive to the derivative
/// of the potential than on the potential itself.  It is therefore
/// absolutely essential that the interpolated derivatives are smooth
/// and continuous.  This function uses simple quadratic interpolation
/// to find f(r).  Since quadric interpolants don't have smooth
/// derivatives, f'(r) is computed using a 4 point finite difference
/// stencil.
///
/// Interpolation is used heavily by the EAM force routine so this
/// function is a potential performance hot spot.  Feel free to
/// reimplement this function (and initInterpolationObject if necessay)
/// with any higher performing implementation of interpolation, as long
/// as the alternate implmentation that has the required smoothness
/// properties.  Cubic splines are one common alternate choice.
///
/// \param [in] table Interpolation table.
/// \param [in] r Point where function value is needed.
/// \param [out] f The interpolated value of f(r).
/// \param [out] df The interpolated value of df(r)/dr.
void interpolateOrig(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{
   const real_t* tt = table->values; // alias
   //real_t tt[] = {0.000039, 0.000038, 0.000037, 0.000036};

   r = (r-table->x0)*(table->invDx) ;
   int ii = (int)floor(r);

   r = r - floor(r);

   real_t g1 = tt[ii+1] - tt[ii-1];
   real_t g2 = tt[ii+2] - tt[ii];

   *f = tt[ii] + 0.5*r*(g1 + r*(tt[ii+1] + tt[ii-1] - 2.0*tt[ii]) );
   *df = 0.5*(g1 + r*(g2-g1))*table->invDx;
}

void interpolateEmpty(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{
   *f = 0.000033;
   *df = -0.000114;
   return;
}

void interpolateFakeTable(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{

  volatile int q = 58;

   real_t tt[] = {0.000039, 0.000038, 0.000037, 0.000036};

   r = (r-table->x0)*(table->invDx) ;
   int ii = 1;

   r = r - floor(r);

   real_t g1 = tt[ii+1] - tt[ii-1];
   real_t g2 = tt[ii+2] - tt[ii];

   *f = tt[ii] + 0.5*r*(g1 + r*(tt[ii+1] + tt[ii-1] - 2.0*tt[ii]) );
   *df = 0.5*(g1 + r*(g2-g1))*table->invDx;
}

void interpolateLocalTable(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{

   //real_t tt[] = {0.000039, 0.000038, 0.000037, 0.000036};
   volatile real_t tt_i_neg_1 = 0.000039;
   volatile real_t tt_i_0 = 0.000038;
   volatile real_t tt_i_1 = 0.000037;
   volatile real_t tt_i_2 = 0.000036;

   r = (r-table->x0)*(table->invDx) ;
   int ii = 1;

   r = r - floor(r);

   real_t g1 = tt_i_1 - tt_i_neg_1;
   real_t g2 = tt_i_2 - tt_i_0;

   *f = tt_i_0 + 0.5*r*(g1 + r*(tt_i_1 + tt_i_neg_1 - 2.0*tt_i_0) );
   *df = 0.5*(g1 + r*(g2-g1))*table->invDx;
}

// this gcc requires movq, clang moves movsd
#ifdef __GNUC__
void interpolateAsmLoad(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{
  const real_t* tt = table->values; // alias
/*
  r = (r-table->x0)*(table->invDx) ;
  int ii = (int)floor(r); */

  int ii;
  asm("movsd  %0, %%xmm0"
     :     /* no output operands */
     :"X"(table->x0));
  asm("movsd %0, %%xmm2" :       :"X"(table->invDx));
  asm("movq %%xmm1, %0" :"=X"(r):                 );
  asm("movl %1, %0"     : "=X"(ii):"X"(100));
  asm("movq  %0, %%xmm2":         :"X"(r));

/*
  // reset r to fractional distance
  r = r - floor(r);*/

  asm("movq %0, %%xmm3"::"X"(r));asm("movq %0, %%xmm4"::"X"(r));asm("movq %%xmm4, %0"::"X"(r));

/*
 real_t g1 = tt[ii+1] - tt[ii-1];
 real_t g2 = tt[ii+2] - tt[ii];
 */

  real_t g1;
  real_t g2;
  asm("movq %1, %0 #added assembly":"=X"(r)        :"X"(tt[ii+1]));
  asm("movsd %0, %%xmm6":        :"X"(tt[ii-1]));
  asm("movq %%xmm6, %0":"=X"(g1):             );
  asm("movsd %0, %%xmm6":        :"X"(tt[ii+2]));
  asm("movsd %1, %0"    :"=X"(g2):"X"(tt[ii  ]));

/*
 *f = tt[ii] + 0.5*r*(g1 + r*(tt[ii+1] + tt[ii-1] - 2.0*tt[ii]) );
 *df = 0.5*(g1 + r*(g2-g1))*table->invDx;
 */

  real_t f_val = 0.000033;
  asm("movsd %0, %%xmm7 #added assembly":     :"X"(tt[ii-1]));
  asm("movsd %0, %%xmm8":     :"X"(tt[ii]));
  asm("movsd %0, %%xmm9":     :"X"(tt[ii+1]));
  asm("movq %0, %%xmm10":     :"X"(r));
  asm("movq %0, %%xmm11":     :"X"(g1));
  asm("movq %0, %%xmm11":     :"X"(g1));
  asm("movsd %0, %%xmm12":     :"X"(table->invDx));
  asm("movq %1, %0":"=X"(*f):"X"(f));

  // to load a usable value
  //*f = 0.000033;
  *df = -0.000114;
  return;
}
#else
void interpolateAsmLoad(InterpolationObject* table, real_t r, real_t* f, real_t* df)
{
   const real_t* tt = table->values; // alias
/*
   r = (r-table->x0)*(table->invDx) ;
   int ii = (int)floor(r); */

   int ii;
   asm("movsd  %0, %%xmm0"
      :     /* no output operands */
      :"X"(table->x0));
   asm("movsd %0, %%xmm2" :       :"X"(table->invDx));
   asm("movsd %%xmm1, %0" :"=X"(r):                 );
   asm("movl %1, %0"     : "=X"(ii):"X"(100));
   asm("movsd  %0, %%xmm2":         :"X"(r));

/*
   // reset r to fractional distance
   r = r - floor(r);*/

   asm("movsd %0, %%xmm3"::"X"(r));asm("movsd %0, %%xmm4"::"X"(r));asm("movsd %%xmm4, %0"::"X"(r));

/*
  real_t g1 = tt[ii+1] - tt[ii-1];
  real_t g2 = tt[ii+2] - tt[ii];
  */

   real_t g1;
   real_t g2;
   asm("movsd %1, %0 #added assembly":"=X"(r)        :"X"(tt[ii+1]));
   asm("movsd %0, %%xmm6":        :"X"(tt[ii-1]));
   asm("movsd %%xmm6, %0":"=X"(g1):             );
   asm("movsd %0, %%xmm6":        :"X"(tt[ii+2]));
   asm("movsd %1, %0"    :"=X"(g2):"X"(tt[ii  ]));

/*
  *f = tt[ii] + 0.5*r*(g1 + r*(tt[ii+1] + tt[ii-1] - 2.0*tt[ii]) );
  *df = 0.5*(g1 + r*(g2-g1))*table->invDx;
  */


   asm("movsd %0, %%xmm7 #added assembly":     :"X"(tt[ii-1]));
   asm("movsd %0, %%xmm8":     :"X"(tt[ii]));
   asm("movsd %0, %%xmm9":     :"X"(tt[ii+1]));
   asm("movsd %0, %%xmm10":     :"X"(r));
   asm("movsd %0, %%xmm11":     :"X"(g1));
   asm("movsd %0, %%xmm11":     :"X"(g1));
   asm("movsd %0, %%xmm12":     :"X"(table->invDx));
   asm("movsd %1, %0":"=X"(*f):"X"(0.000033));

   // to load a usable value
   //*f = 0.000033;
   *df = -0.000114;
}
#endif
