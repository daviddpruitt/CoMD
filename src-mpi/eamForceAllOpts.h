/*
 * eamForce kernel with neighborlists and seperate loops
 * selectable through defines, makes code harder to read though
 * defaults to the standard implementation
 */
int eamForce(SimFlat* s)
{
   EamPotential* pot = (EamPotential*) s->pot;
   assert(pot);

   // set up halo exchange and internal storage on first call to forces.
   if (pot->forceExchange == NULL)
   {
      int maxTotalAtoms = MAXATOMS*s->boxes->nTotalBoxes;
      pot->dfEmbed = comdMalloc(maxTotalAtoms*sizeof(real_t));
      pot->rhobar  = comdMalloc(maxTotalAtoms*sizeof(real_t));
      pot->forceExchange = initForceHaloExchange(s->domain, s->boxes);
      pot->forceExchangeData = comdMalloc(sizeof(ForceExchangeData));
      pot->forceExchangeData->dfEmbed = pot->dfEmbed;
      pot->forceExchangeData->boxes = s->boxes;
   }

   real_t rCut2 = pot->cutoff*pot->cutoff;

   // zero forces / energy / rho /rhoprime
   real_t etot = 0.0;
   memset(s->atoms->f,  0, s->boxes->nTotalBoxes*MAXATOMS*sizeof(real3));
   memset(s->atoms->U,  0, s->boxes->nTotalBoxes*MAXATOMS*sizeof(real_t));
   memset(pot->dfEmbed, 0, s->boxes->nTotalBoxes*MAXATOMS*sizeof(real_t));
   memset(pot->rhobar,  0, s->boxes->nTotalBoxes*MAXATOMS*sizeof(real_t));

   // prefetch pots
   prefetchPots(pot);

   int nbrBoxes[27];
   // loop over local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int nIBox = s->boxes->nAtoms[iBox];
      int nNbrBoxes = getNeighborBoxes(s->boxes, iBox, nbrBoxes);
      // loop over neighbor boxes of iBox (some may be halo boxes)
      for (int jTmp=0; jTmp<nNbrBoxes; jTmp++)
      {
         int jBox = nbrBoxes[jTmp];
         if (jBox < iBox ) continue;

	       prefetchBoxes(iBox, jBox, s);
	       frivStartTimer(computeForce1Timer);

         int nJBox = s->boxes->nAtoms[jBox];
         // loop over atoms in iBox
         for (int iOff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,iOff++)
         {
#ifdef NEIGHBOR_LIST
            // set # of neighbors in neighbor box
            int numNeighbors = s->atoms->numNeighbors[iOff];
            real_t dr[MAXATOMS];
            real_t r2[MAXATOMS];// {0};

            // loop over atoms in jBox
            for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
            {
              real_t drx = s->atoms->r[iOff][0]-s->atoms->r[jOff][0];
              r2[ij]=drx*drx;
            }

            // compute vector differences
            for (int k=1; k<3; k++)
            {
              // loop over atoms in jBox
              for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
              {
                real_t drx = s->atoms->r[iOff][k]-s->atoms->r[jOff][k];
                r2[ij]+=drx*drx;
              }
            }

            /* real_t dr[MAXATOMS][3]; */
            /* real_t r2[MAXATOMS] = {0}; */

            /* // loop over atoms in jBox */
            /* 	for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++) */
            /* 	{ */
            /*   // compute vector differences */
            /*   for (int k=0; k<3; k++) */
            /*   { */
	    /*           dr[k][ij] = s->atoms->r[iOff][k]-s->atoms->r[jOff][k]; */
	    /*           r2[ij]+=dr[k][ij]*dr[k][ij]; */
            /*   } */
            /* } */

            // loop over atoms in jBox
            for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
            {
               if ( (iBox==jBox) &&(ij <= ii) ) continue;

               if(r2[ij]<rCut2) {
                 // if we're within cutoff add to neighbors
                 s->atoms->neighbors[iOff].neighborList[numNeighbors] = jOff;
                 //s->atoms->neighbors[iOff].distance[numNeighbors] = sqrt(r2);
                 if (jBox < s->boxes->nLocalBoxes)
                   s->atoms->neighbors[iOff].neighborScale[numNeighbors] = 2;
                 else
                   s->atoms->neighbors[iOff].neighborScale[numNeighbors] = 1;

                 numNeighbors++;
               }
            } // loop over atoms in jBox
            s->atoms->numNeighbors[iOff] = numNeighbors;
            // if (max_neighbors < numNeighbors)
            //   max_neighbors = numNeighbors;
#else
            for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
            {
               if ( (iBox==jBox) &&(ij <= ii) ) continue;

               double r2 = 0.0;
               real3 dr;
               for (int k=0; k<3; k++)
               {
                  dr[k]=s->atoms->r[iOff][k]-s->atoms->r[jOff][k];
                  r2+=dr[k]*dr[k];
               }
               if(r2>=rCut2) continue;
               real_t r = sqrt(r2);

               real_t phiTmp, dPhi, rhoTmp, dRho;
               interpolate(pot->phi, r, &phiTmp, &dPhi);
               interpolate(pot->rho, r, &rhoTmp, &dRho);

               for (int k=0; k<3; k++)
               {
                  s->atoms->f[iOff][k] -= dPhi*dr[k]/r;
                  s->atoms->f[jOff][k] += dPhi*dr[k]/r;
               }

               // update energy terms
               // calculate energy contribution based on whether
               // the neighbor box is local or remote
               if (jBox < s->boxes->nLocalBoxes)
                  etot += phiTmp;
               else
                  etot += 0.5*phiTmp;

               s->atoms->U[iOff] += 0.5*phiTmp;
               s->atoms->U[jOff] += 0.5*phiTmp;

               // accumulate rhobar for each atom
               pot->rhobar[iOff] += rhoTmp;
               pot->rhobar[jOff] += rhoTmp;

            } // loop over atoms in jBox
#endif
         } // loop over atoms in iBox
	 frivStopTimer(computeForce1Timer);

      } // loop over neighbor boxes
   } // loop over local boxes

#ifdef NEIGHBOR_LIST
   // loop over neighborLists
   // loop over local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int nIBox = s->boxes->nAtoms[iBox];
      // loop over atoms in ibox
       for (int iOff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,iOff++)
       {
          // temp structures
          Atoms *restrict atomsI = s->atoms;
          Atoms *restrict atomsJ = s->atoms;
          int numNeighbors = atomsI->numNeighbors[iOff];

          real3 dr[MAX_NEIGHBORS];
          real_t phiTmp[MAX_NEIGHBORS], dPhi[MAX_NEIGHBORS];
          real_t rhoTmp[MAX_NEIGHBORS], dRho[MAX_NEIGHBORS];
          real_t r[MAX_NEIGHBORS];
          real_t r2[MAX_NEIGHBORS];

          for (int i = 0; i < numNeighbors; ++i) {
             int jOff = atomsI->neighbors[iOff].neighborList[i];
# ifndef SEP_LOOPS
             r2[i] = 0.0;
             // get distances
             for (int k=0; k<3; k++)
             {
                dr[i][k]=atomsI->r[iOff][k]-atomsJ->r[jOff][k];
                r2[i]+=dr[i][k]*dr[i][k];
             }
# endif
# ifdef SEP_LOOPS
             // get distances
             for (int k=0; k<3; k++)
             {
                dr[i][k]=atomsI->r[iOff][k]-atomsJ->r[jOff][k];
             }
          }
          for (int i = 0; i < numNeighbors; ++i) {
             r2[i] = 0.0;
             for (int k=0; k<3; k++)
             {
                r2[i]+=dr[i][k]*dr[i][k];
             }
          }

          for (int i = 0; i < numNeighbors; ++i) {
             int jOff = atomsI->neighbors[iOff].neighborList[i];
# endif
             r[i] = sqrt(r2[i]);//atomsI->neighbors[iOff].distance[i];
             interpolate(pot->phi, r[i], &phiTmp[i], &dPhi[i]);
             interpolate(pot->rho, r[i], &rhoTmp[i], &dRho[i]);
# ifdef SEP_LOOPS
          }
          for (int i = 0; i < numNeighbors; ++i) {
             int jOff = atomsI->neighbors[iOff].neighborList[i];
# endif
             for (int k=0; k<3; k++)
             {
                atomsI->f[iOff][k] -= dPhi[i]*dr[i][k]/r[i];
                atomsJ->f[jOff][k] += dPhi[i]*dr[i][k]/r[i];
             }

             // update energy terms
             // calculate energy contribution based on whether
             // the neighbor box is local or remote
             etot += (0.5*atomsI->neighbors[iOff].neighborScale[i])*phiTmp[i];

             atomsI->U[iOff] += 0.5*phiTmp[i];
             atomsJ->U[jOff] += 0.5*phiTmp[i];

             // accumulate rhobar for each atom
             pot->rhobar[iOff] += rhoTmp[i];
             pot->rhobar[jOff] += rhoTmp[i];
          } // loop over neighbors
          // done with the atoms in this box, reset it's counters
          //s->atoms->numNeighbors[iOff] = 0;
       } // loop over atoms in iBox
	 frivStopTimer(computeForce1Timer);
   } // loop over local boxes
#endif

   // Compute Embedding Energy
   // loop over all local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int iOff;
      int nIBox =  s->boxes->nAtoms[iBox];

      // loop over atoms in iBox
      frivStartTimer(computeForce2Timer);
      for (int iOff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,iOff++)
      {
         real_t fEmbed, dfEmbed;
         interpolate(pot->f, pot->rhobar[iOff], &fEmbed, &dfEmbed);
         pot->dfEmbed[iOff] = dfEmbed; // save derivative for halo exchange
         etot += fEmbed;
         s->atoms->U[iOff] += fEmbed;
      }
      frivStopTimer(computeForce2Timer);
   }

   // exchange derivative of the embedding energy with repsect to rhobar
   startTimer(eamHaloTimer);
   haloExchange(pot->forceExchange, pot->forceExchangeData);
   stopTimer(eamHaloTimer);

   // third pass
#ifdef NEIGHBOR_LIST
   // loop over local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int nIBox = s->boxes->nAtoms[iBox];
      // loop over atoms in ibox
       for (int iOff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,iOff++)
       {
          // loop over neighborLists
          // temp structures
          Atoms *restrict atomsI = s->atoms;
          Atoms *restrict atomsJ = s->atoms;
          int numNeighbors = atomsI->numNeighbors[iOff];
          real_t r[MAX_NEIGHBORS];
          real_t rhoTmp[MAX_NEIGHBORS], dRho[MAX_NEIGHBORS];
          real3 dr[MAX_NEIGHBORS];
          real_t r2[MAX_NEIGHBORS];

          // loop over neighbors
          for (int i = 0; i < numNeighbors; ++i) {
             int jOff = atomsI->neighbors[iOff].neighborList[i];
# ifndef SEP_LOOPS
             r2[i] = 0.0;
             // get distances
             for (int k=0; k<3; k++)
             {
                dr[i][k]=atomsI->r[iOff][k]-atomsJ->r[jOff][k];
                r2[i]+=dr[i][k]*dr[i][k];
             }
# endif
# ifdef SEP_LOOPS
             // get distances
             for (int k=0; k<3; k++)
             {
                dr[i][k]=atomsI->r[iOff][k]-atomsJ->r[jOff][k];
             }
          }
          for (int i = 0; i < numNeighbors; ++i) {
             r2[i] = 0.0;
             for (int k=0; k<3; k++)
             {
                r2[i]+=dr[i][k]*dr[i][k];
             }
          }
          for (int i = 0; i < numNeighbors; ++i) {
# endif
             r[i] = sqrt(r2[i]);//atomsI->neighbors[iOff].distance[i];
             interpolate(pot->rho, r[i], &rhoTmp[i], &dRho[i]);
# ifdef SEP_LOOPS
          }
          for (int i = 0; i < numNeighbors; ++i) {
             int jOff = atomsI->neighbors[iOff].neighborList[i];
# endif
             for (int k=0; k<3; k++)
             {
                atomsI->f[iOff][k] -= (pot->dfEmbed[iOff]+pot->dfEmbed[jOff])*dRho[i]*dr[i][k]/r[i];
                atomsJ->f[jOff][k] += (pot->dfEmbed[iOff]+pot->dfEmbed[jOff])*dRho[i]*dr[i][k]/r[i];
             }
          } // loop over neighbors
          // done with the atoms in this box, reset it's counters
          s->atoms->numNeighbors[iOff] = 0;
       } // loop over atoms in iBox
	 frivStopTimer(computeForce1Timer);
   } // loop over local boxes
#else
   // loop over local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int nIBox =  s->boxes->nAtoms[iBox];
      int nNbrBoxes = getNeighborBoxes(s->boxes, iBox, nbrBoxes);
      // loop over neighbor boxes of iBox (some may be halo boxes)
      for (int jTmp=0; jTmp<nNbrBoxes; jTmp++)
      {
         int jBox = nbrBoxes[jTmp];
         if(jBox < iBox) continue;

	       prefetchBoxes(iBox, jBox, s);
	       frivStartTimer(computeForce3Timer);

         int nJBox = s->boxes->nAtoms[jBox];
         // loop over atoms in iBox
         for (int iOff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,iOff++)
         {
            // loop over atoms in jBox
            for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
            {
               if ((iBox==jBox) && (ij <= ii))  continue;

               double r2 = 0.0;
               real3 dr;
               for (int k=0; k<3; k++)
               {
                  dr[k]=s->atoms->r[iOff][k]-s->atoms->r[jOff][k];
                  r2+=dr[k]*dr[k];
               }
               if(r2>=rCut2) continue;

               real_t r = sqrt(r2);

               real_t rhoTmp, dRho;
               interpolate(pot->rho, r, &rhoTmp, &dRho);

               real3 *restrict atomsI = &s->atoms->f[iOff];
               real3 *restrict atomsJ = &s->atoms->f[jOff];

               for (int k=0; k<3; k++)
               {
                  s->atoms->f[iOff][k] -= (pot->dfEmbed[iOff]+pot->dfEmbed[jOff])*dRho*dr[k]/r;
                  s->atoms->f[jOff][k] += (pot->dfEmbed[iOff]+pot->dfEmbed[jOff])*dRho*dr[k]/r;
               }

            } // loop over atoms in jBox
         } // loop over atoms in iBox
	 frivStopTimer(computeForce3Timer);
      } // loop over neighbor boxes
   } // loop over local boxes
#endif
   s->ePotential = (real_t) etot;

   return 0;
}
