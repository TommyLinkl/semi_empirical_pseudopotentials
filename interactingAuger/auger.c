#include "fd.h"

/*****************************************************************************/

void calculate_auger(double *Cbs,double *Ebs,double *vkijb,double *vkcab,double *evalbe,double *evalai,double *sigeai,lng_st ist,par_st par)
{
  FILE *pf; 
  long a, i, b, c, j, k, jms, ibs, *list, *lista, *listi;  
  long iBiexc, iExc1, iExc2, thl, th2l, thl2, numConsWindows = 10;
  double *w2h, *w2e, tmp;
  double *eRate, *hRate, *biexcEnergy;
  double sum, pF, this_lifetime_ps;
  long *eCount, *hCount;

  // Write time this function began
  writeFunctionStartTime("Program is now calculating the auger recombination rates:", stdout);

  // Allocate memory used within this function
  if ((w2h = (double*)calloc(ist.numBandEdgeElectrons*ist.itot*ist.numBiexcitons,sizeof(double)))==NULL) nerror("w2h");
  if ((w2e = (double*)calloc(ist.atot*ist.numBandEdgeHoles*ist.numBiexcitons,sizeof(double)))==NULL) nerror("w2e");
  if ((eRate = (double*)calloc(ist.numBiexcitons*numConsWindows,sizeof(double)))==NULL) nerror("eRate");
  if ((hRate = (double*)calloc(ist.numBiexcitons*numConsWindows,sizeof(double)))==NULL) nerror("hRate");
  if ((eCount = (long*)calloc(ist.numBiexcitons*numConsWindows,sizeof(double)))==NULL) nerror("eCount");
  if ((hCount = (long*)calloc(ist.numBiexcitons*numConsWindows,sizeof(double)))==NULL) nerror("hCount");
  if ((biexcEnergy = (double*)calloc(ist.numBiexcitons,sizeof(double)))==NULL) nerror("biexcEnergy");

  if ((list = (long*)calloc(ist.numBandEdgeElectrons*ist.numBandEdgeHoles,sizeof(long)))==NULL) nerror("list");
  for (ibs = 0, a = 0; a < ist.numBandEdgeElectrons; a++)
    for (i = 0; i < ist.numBandEdgeHoles; i++, ibs++) list[a*ist.numBandEdgeHoles+i] = ibs;

  if ((listi = (long*)calloc(ist.itot,sizeof(long)))==NULL) nerror("list i");
  for (i = 0, jms = 0; jms < ist.ms; jms++) 
    if ((evalai[jms] >= par.Eimin) && (evalai[jms] <= par.Eimax) && (sigeai[jms] < par.deps)) {listi[i] = jms; i++;}
  
  if ((lista = (long*)calloc(ist.atot,sizeof(long)))==NULL) nerror("list a");
  for (a = 0, jms = 0; jms < ist.ms; jms++) 
    if ((evalai[jms] >= par.Eamin) && (evalai[jms] <= par.Eamax) && (sigeai[jms] < par.deps)) {lista[a] = jms; a++;}
  
  // Calculate w2h 
  thl = ist.numBandEdgeHoles * ist.numBandEdgeElectrons;
  th2l = thl * ist.numBandEdgeHoles;
  iBiexc = 0;
  for (iExc1 = 0; iExc1 < ist.numExcitons; iExc1++) {
    for (iExc2 = 0; iExc2 < ist.numExcitons; iExc2++) { 
      biexcEnergy[iBiexc] = Ebs[iExc1] + Ebs[iExc2];
      if (biexcEnergy[iBiexc] < par.maxInitE+EPS) {
        for (a = 0; a < ist.numBandEdgeElectrons; a++) { 
          for (i = 0; i < ist.itot; i++) {
            for (sum = 0.0, k = 0; k < ist.numBandEdgeHoles; k++) {
              for (j = 0; j < ist.numBandEdgeHoles; j++) {
                for (b = 0; b < ist.numBandEdgeElectrons; b++) {  
                  sum += vkijb[i*th2l + j*thl + k*ist.numBandEdgeElectrons + b]  
                    * Cbs[iExc1*ist.msbs2 + list[b*ist.numBandEdgeHoles+j]] * Cbs[iExc2*ist.msbs2 + list[a*ist.numBandEdgeHoles+k]];
                 }
               }
             }
           w2h[iBiexc*ist.numBandEdgeElectrons*ist.itot + a*ist.itot + i] = sum*sum;
           }
         }
         iBiexc++;
       }
     }
   }

  // Calculate w2e
  thl = ist.numBandEdgeHoles * ist.numBandEdgeElectrons;
  thl2 = thl * ist.numBandEdgeElectrons;
  iBiexc = 0;
  for (iExc1 = 0; iExc1 < ist.numExcitons; iExc1++) {
    for (iExc2 = 0; iExc2 < ist.numExcitons; iExc2++) {
      biexcEnergy[iBiexc] = Ebs[iExc1] + Ebs[iExc2];
      if (biexcEnergy[iBiexc] < par.maxInitE+EPS) {
        for (a = 0; a < ist.atot; a++) {
          for (i = 0; i < ist.numBandEdgeHoles; i++) {
            for (sum = 0.0, b = 0; b < ist.numBandEdgeElectrons; b++) {
              for (k = 0; k < ist.numBandEdgeHoles; k++) {
                for (c = 0; c < ist.numBandEdgeElectrons; c++) {
                  sum += vkcab[a*thl2 + b*thl + k*ist.numBandEdgeElectrons + c] 
                    * Cbs[iExc1*ist.msbs2 + list[b*ist.numBandEdgeHoles+i]] * Cbs[iExc2*ist.msbs2 + list[c*ist.numBandEdgeHoles+k]];
                 }
               }
             }     
           w2e[iBiexc*ist.atot*ist.numBandEdgeHoles + a*ist.numBandEdgeHoles + i] = sum*sum;
           }
         }
         iBiexc++;
       }
     }
   }  

  // Create multiple energy conservations windows so do not have to repeat calculation
  long iDeltaE; 
  double deltaE[numConsWindows]; // energy conservation windows -> equally spaced from par.DeltaE to 0.1*par.DeltaE
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) deltaE[iDeltaE] = par.DeltaE * (1.0 - 0.1*(double)(iDeltaE));

  // Calculate k_x+
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.numBiexcitons; iBiexc++) {
      for (a = 0; a < ist.numBandEdgeElectrons; a++) {
        for (i = 0; i < ist.itot; i++) {
		  if (fabs(evalbe[ist.lumoIndex + a] - evalai[list[i]] - biexcEnergy[iBiexc]) <= deltaE[iDeltaE]) { // energy conservation check
            hRate[iBiexc + iDeltaE*ist.numBiexcitons] += w2h[iBiexc*ist.numBandEdgeElectrons*ist.itot + a*ist.itot + i];
            hCount[iBiexc + iDeltaE*ist.numBiexcitons]++;
          }
		}
      }
    }  
  }
 
  // Calculate k_x-
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.numBiexcitons; iBiexc++) {  
      for (a = 0; a < ist.atot; a++) {
        for (i = 0; i < ist.numBandEdgeHoles; i++) {
          if (fabs(evalai[lista[a]] - evalbe[i] - biexcEnergy[iBiexc]) <= deltaE[iDeltaE]) { // energy conservation check
            eRate[iBiexc + iDeltaE*ist.numBiexcitons] += w2e[iBiexc*ist.atot*ist.numBandEdgeHoles + a*ist.numBandEdgeHoles + i];
            eCount[iBiexc + iDeltaE*ist.numBiexcitons]++; 
          }
        }
	  }
    }
  }

  // rates scaled to include state degeneracy and finite width of the energy conservation delta function
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.numBiexcitons; iBiexc++ ) {  
      hRate[iBiexc + iDeltaE*ist.numBiexcitons] *= TWOPI/(2.0*deltaE[iDeltaE]); 
      eRate[iBiexc + iDeltaE*ist.numBiexcitons] *= TWOPI/(2.0*deltaE[iDeltaE]); 
    }
  }

  // prints the auger rate for each given initial biexcitonic sate 
  // prints out: iExc1, iExc2, biexcEnergy, numHotElecs, numHotHoles, 
  // k_x- (au & ps), k_x+ (au & ps), k_AR (au & ps), lifetime(ps)
  // TODO: add loop over numConsWindows -> different file names for each
  pf = fopen("auger.dat" , "w");
  fprintf (pf,"# iExc1, iExc2, biexcEnergy, numHotElecs, numHotHoles, k_x- (au & ps), k_x+ (au & ps), k_AR (au & ps), lifetime(ps)\n");      
  iBiexc = 0;
  for (iExc1 = 0; iExc1 < ist.numExcitons; iExc1++) {
    for (iExc2 = 0; iExc2 < ist.numExcitons; iExc2++) {
      biexcEnergy[iBiexc] = Ebs[iExc1] + Ebs[iExc2];
      if (biexcEnergy[iBiexc] < par.maxInitE+EPS) {
        fprintf (pf,"%ld %ld %.10f %ld %ld %g %g %g %g %g %g %g\n", iExc1, iExc2, biexcEnergy[iBiexc], 
          eCount[iBiexc], hCount[iBiexc], eRate[iBiexc], eRate[iBiexc]*AUTOPS, hRate[iBiexc], hRate[iBiexc]*AUTOPS, 
          eRate[iBiexc]+hRate[iBiexc], (eRate[iBiexc]+hRate[iBiexc])*AUTOPS, 1.0/((eRate[iBiexc]+hRate[iBiexc])*AUTOPS));      
        iBiexc++;
      }
    }
  }
  fclose(pf);

  // Perform Boltzmann weighted average over the biexciton initial states
  pf = fopen("window_augerLifetime.dat", "w");
  fprintf(pf, "# Energy_Window (a.u.)       Auger_Lifetime (ps) \n");
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    if (iDeltaE) writeSeparation(stdout);
    printf("Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    this_lifetime_ps = calc_boltzmann_weighted_rates(biexcEnergy, &(eRate[iDeltaE*ist.numBiexcitons]), 
      &(hRate[iDeltaE*ist.numBiexcitons]), par.temp, ist.numBiexcitons);
    fprintf(pf, "%.6f      %.12f \n", deltaE[iDeltaE], this_lifetime_ps);
  }
  fclose(pf); 

  // Free dynamically allocated memory 
  free(w2h); free(w2e); 
  free(lista); free(listi); free(list);
  free(eRate); free(hRate); free(eCount); free(hCount); 
  free(biexcEnergy);

  return;
}

/*****************************************************************************/
// calculates and prints the boltzmann weighted AR, k_x-, k_x+ rates and 
// the corresponding lifetimes and boltzmann populations

double calc_boltzmann_weighted_rates(double *energies, double *eRate, double *hRate, double temp, int numStates) {
  FILE *pf;
  long i;
  double pF, bwAugerRate, bwElecRate, bwHoleRate;
  double stateProb, e0 = energies[0];  // will be energy of the lowest state
  double beta = AUTOEV/(KB*(temp+EPS));  

  bwAugerRate = bwElecRate = bwHoleRate = 0.0;
  pf = fopen("boltzmannStats.dat", "w");
  pF = calc_partition_function(energies, temp, numStates);

  for (i = 0; i < numStates; i++) if (energies[i] < e0) e0 = energies[i]; // in case energies isn't ordered
  for (i = 0; i < numStates; i++) {
    stateProb = exp(-beta*(energies[i]-e0));
    bwElecRate += eRate[i]*stateProb;
    bwHoleRate += hRate[i]*stateProb;
    fprintf(pf, "%.10f %.8f %.12f %.12f %.12f\n", energies[i], stateProb/pF, 
      1.0/eRate[i]/AUTOPS, 1.0/hRate[i]/AUTOPS, 1.0/(eRate[i]+hRate[i])/AUTOPS);
  }
  fclose(pf);

  bwElecRate /= pF;
  bwHoleRate /= pF;
  bwAugerRate = bwElecRate + bwHoleRate;

  printf("Lowest energy biexciton = %.12f\n", e0);
  printf("Partition function = %.12f\n", pF);
  printf("Boltzmann weighted k_x- = %.12f %.6f\n", bwElecRate, bwElecRate*AUTOPS);
  printf("Boltzmann weighted k_x+ = %.12f %.6f\n", bwHoleRate, bwHoleRate*AUTOPS);
  printf("Boltzmann weighted kAR  = %.12f %.6f\n", bwAugerRate, bwAugerRate*AUTOPS);
  printf("Boltzmann weighted x- lifetime = %.6f a.u. %.12f ps\n", 1.0/bwElecRate, 1.0/(bwElecRate*AUTOPS));
  printf("Boltzmann weighted x+ lifetime = %.6f a.u. %.12f ps\n", 1.0/bwHoleRate, 1.0/(bwHoleRate*AUTOPS));
  printf("Boltzmann weighted Auger lifetime = %.6f a.u. %.12f ps\n", 1.0/(bwAugerRate), 1.0/(bwAugerRate*AUTOPS));

  return (1.0/(bwAugerRate*AUTOPS));
}

/*****************************************************************************/
// returns the partition function given an array of energies, a temperature,
// and the number of states (i.e., length of the energies array)
// no state degeneracy is included here

double calc_partition_function(double *energies, double temp, int numStates) {
  long i;
  double pF = 0.0;
  double e0 = energies[0];  // will be energy of the lowest state
  double beta = AUTOEV/(KB*(temp+EPS)); 

  for (i = 0; i < numStates; i++) if (energies[i] < e0) e0 = energies[i]; // in case energies isn't ordered
  for (i = 0; i < numStates; i++) pF += exp(-beta*(energies[i]-e0));

  return pF;
}

/*****************************************************************************/
