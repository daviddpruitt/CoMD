/// Interpolate a table to determine f(r) and its derivative f'(r).for a vector of distances
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
/// \param [in] rvIn Points r where function value is needed.
/// \param [out] fv The interpolated values of f(r).
/// \param [out] dfv The interpolated values of df(r)/dr.
/// \param [in] vlen Number of points to process
#define imin(_a, _b) ( ((_a) <= (_b)) ? (_a) : (_b) ) /* integer min */
#define VLEN 64			// strip mine up to 32 at a time
#pragma omp declare simd
void interpolatev(InterpolationObject *restrict table, real_t *restrict rvIn, real_t *restrict fv, real_t *restrict dfv, int vlen)
{
#ifdef FAKE_V
  static int once = 0;
  for (int i = 0; i < vlen; i++) 
    interpolate(table, rvIn[i], &fv[i], &dfv[i]);
  return;
#else // ndef FAKE_V

  const real_t *restrict tt = table->values; /* alias */
  const real_t x0 = table->x0, invDx = table->invDx; /* alias */
  real_t numSamples = (float)table->n;
  int chunkSize = 0; 		                      // chunksize computed inside loop
  

  for (int chunkBegin = 0; chunkBegin < vlen; chunkBegin += chunkSize) { /* chunk */
    real_t g1v[VLEN], g2v[VLEN], samplv[4][VLEN], rv[VLEN]; // temporaries for vector ops
    int iiv[VLEN];		                            // indices into table
    const int chunkLimit = imin(chunkBegin+VLEN, vlen);   // determine chunk range
    chunkSize = chunkLimit - chunkBegin;              //  # elements in chunk
#pragma vector aligned
#pragma omp simd
    for (int i = 0, ir = chunkBegin; i < chunkSize; i++, ir++)               // normalize to table indices (Vectorizable)
      rv[i] = fmin(fmax((rvIn[ir]-x0)*invDx, 0.0), (real_t)numSamples);
#pragma vector aligned
    for (int i = 0; i < chunkSize; i++) {             // compute indices for table slices (Vectorizable)
      long floorR = rv[i];
      iiv[i] = floorR - 1;		      // index -1 is ok because values table is index-biased + 1
      rv[i] = rv[i] - floorR;			      // keep fractional part
      //      { printf("interpolatev: inR[0]=%f, rv[0]=%f, iiv[i]=%d\n", rvIn[0], rv[0], iiv[0]); fflush(stdout); exit(1); }
    }
    for (int i = 0; i < chunkSize; i++) {             // scatter table slices
      for (int j = 0, ii = iiv[i]; j < 4; j++, ii++) 
	samplv[j][i] = tt[ii];
    }
/* #pragma vector aligned */
/* #pragma omp simd */
/*     for (int i = 0; i < chunkSize; i++) {             // reset r to fractional distance (Vectorizable) */
/*       rv[i] = rv[i] - (int)(rv[i]); */
/*     } */
#pragma vector aligned
#pragma omp simd
    for (int i = 0; i < chunkSize; i++) {             // compute g1v & g2v (Vectorizable)
      //g1   = tt[ii+1]  -    tt[ii-1];
      g1v[i] = samplv[2][i] - samplv[0][i];
      // g2 =  tt[ii+2]     - tt[ii];
      g2v[i] = samplv[3][i] - samplv[1][i];
    }
#pragma vector aligned
#pragma omp simd
    for (int i = 0, ir = chunkBegin; i < chunkSize; i++, ir++) {             // compute f & df
      //*f  = tt[ii]       + 0.5*r    *(g1     + r    *(tt[ii+1]     + tt[ii-1]     - 2.0*tt[ii]      ) );
      fv[ir] = samplv[1][i] + 0.5*rv[i]*(g1v[i] + rv[i]*(samplv[2][i] + samplv[0][i] - 2.0*samplv[1][i]) );

      // *df = 0.5*(g1     + r    *(g2     - g1))    *table->invDx;
      dfv[ir] = 0.5*(g1v[i] + rv[i]*(g2v[i] - g1v[i]))*invDx;
      //      { printf("interpolatev: inR[0]=%f, rv[0]=%f, iiv[i]=%d\n", rvIn[0], rv[0], iiv[0]); fflush(stdout); exit(1); }
    }
  } /* chunk */
#endif // ndef FAKE_V
}
