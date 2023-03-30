#include "fd.h"

/*****************************************************************************/

void calculate_auger(double *vkijb,double *vkcab,double *evalbe,double *evalai,double *sigeai,lng_st ist,par_st par)
{
  FILE *pf; 
  long a, i, b, c, j, k, jms, *lista, *listi;  
  long iBiexc, thl, th2l, thl2, numConsWindows = 10;
  double tmp, *eRate, *hRate, *biexcEnergy;
  long *eCount, *hCount;

  if ((eRate = (double *) calloc(ist.numBiexcitons*numConsWindows, sizeof(double)))==NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.numBiexcitons*numConsWindows, sizeof(double)))==NULL) nerror("hRate");
  if ((eCount = (long *) calloc(ist.numBiexcitons*numConsWindows, sizeof(double)))==NULL) nerror("eCount");
  if ((hCount = (long *) calloc(ist.numBiexcitons*numConsWindows, sizeof(double)))==NULL) nerror("hCount");
  if ((biexcEnergy = (double *) calloc(ist.numBiexcitons, sizeof(double)))==NULL) nerror("biexcEnergy");

  if ((listi = (long*)calloc(ist.itot,sizeof(long)))==NULL) nerror("list i");
  for (i = 0, jms = 0; jms < ist.ms; jms++) 
    if ((evalai[jms] >= par.Eimin) && (evalai[jms] <= par.Eimax) && (sigeai[jms] < par.deps)) {listi[i] = jms; i++;}
  
  if ((lista = (long*)calloc(ist.atot,sizeof(long)))==NULL) nerror("list a");
  for (a = 0, jms = 0; jms < ist.ms; jms++) 
    if ((evalai[jms] >= par.Eamin) && (evalai[jms] <= par.Eamax) && (sigeai[jms] < par.deps)) {lista[a] = jms; a++;}
  
  // Create multiple energy conservation windows so do not have to repeat calculation
  long iDeltaE;
  double deltaE[numConsWindows]; // energy conservation windowns -> equally spaced from par.DeltaE to 0.1*par.DeltaE
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) deltaE[iDeltaE] = par.DeltaE * (1.0 - 0.1*(double)(iDeltaE));

  // Calculate k_x+ and k_x-
  thl = ist.numBandEdgeHoles * ist.numBandEdgeElectrons;
  th2l = thl * ist.numBandEdgeHoles;
  thl2 = thl * ist.numBandEdgeElectrons;
  iBiexc = 0;
  // loop over all energy conservation windows
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    iBiexc = 0;
    // loop over all initial biexcitonic states (b,j,k,c)
    for (b = 0; b < ist.numBandEdgeElectrons; b++) { 
      for (j = 0; j < ist.numBandEdgeHoles; j++) {
        for (k = 0; k < ist.numBandEdgeHoles; k++) {
          for (c = 0; c < ist.numBandEdgeElectrons; c++) {  
            biexcEnergy[iBiexc] = evalbe[ist.lumoIndex+b]+evalbe[ist.lumoIndex+c]-evalbe[j]-evalbe[k];
            if (biexcEnergy[iBiexc] < (par.maxInitE+EPS)) { // checks to see if initial biexciton is in allowed energy range
              // Calculate k_x+
              // loop over all final states (delta_ac so just over i)
              for (i = 0; i < ist.itot; i++) { 
                // check to see if the final state is in energy window for the given initial biexciton energy
                if (fabs(evalbe[ist.lumoIndex+c] - evalai[listi[i]] - biexcEnergy[iBiexc]) <= deltaE[iDeltaE]) {
                  tmp = vkijb[i*th2l + j*thl + k*ist.numBandEdgeElectrons + b];
                  hRate[iBiexc + iDeltaE*ist.numBiexcitons] += tmp*tmp;
                  hCount[iBiexc + iDeltaE*ist.numBiexcitons]++;
                }
              }
              // Calculate k_x-        
              // loop over all final states (delta_ij so just over a)
              for (a = 0; a < ist.atot; a++) { 
                // check to see if state a is in energy window for the given initial biexciton energy
                if (fabs(evalai[lista[a]] - evalbe[j] - biexcEnergy[iBiexc]) <= deltaE[iDeltaE]) {
                  tmp = vkcab[a*thl2 + b*thl + k*ist.numBandEdgeElectrons + c];
                  eRate[iBiexc + iDeltaE*ist.numBiexcitons] += tmp*tmp;
                  eCount[iBiexc + iDeltaE*ist.numBiexcitons]++;
                 }
              }     
              iBiexc++; // increments biexciton index energy range of biexciton was satisfied
            }
          }
        }
      }
    }
  }
  // rates scaled to include state degeneracy and finite width of the energy conservation delta function
  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.numBiexcitons; iBiexc++) {
      hRate[iBiexc + iDeltaE*ist.numBiexcitons] *= TWOPI/(2.0*deltaE[iDeltaE]);
      eRate[iBiexc + iDeltaE*ist.numBiexcitons] *= TWOPI/(2.0*deltaE[iDeltaE]);
    }
  }

  // prints the auger rate for each given initial biexcitonic sate - will need to boltzmann on top of this
  // prints out: ie1, ih1, ie2, ih2, biexcEnergy, numHotElecs, numHotHoles, k_x- (au & ps), k_x+ (au & ps), k_AR (au & ps), lifetime(ps)
  // TODO: add loop over numConsWindowns -> different file names for each. Currently just prints AR rates for par.DeltaE choice
  pf = fopen("auger.dat" , "w");
  iBiexc = 0;
  for (b = 0; b < ist.numBandEdgeElectrons; b++) {
    for (j = 0; j < ist.numBandEdgeHoles; j++) { // will be delta function on j
      for (k = 0; k < ist.numBandEdgeHoles; k++) {
        for (c = 0; c < ist.numBandEdgeElectrons; c++) {  
          biexcEnergy[iBiexc] = evalbe[ist.lumoIndex+b]+evalbe[ist.lumoIndex+c]-evalbe[j]-evalbe[k];
          if (biexcEnergy[iBiexc] < (par.maxInitE+EPS)) { // checks to see if initial biexciton is in allowed energy range
            fprintf (pf,"%ld %ld %ld %ld %.10f %ld %ld %g %g %g %g %g %g %g\n", 
              ist.lumoIndex+b, j, ist.lumoIndex+c, k, 
              biexcEnergy[iBiexc], eCount[iBiexc], hCount[iBiexc], 
              eRate[iBiexc], eRate[iBiexc]*AUTOPS, hRate[iBiexc], hRate[iBiexc]*AUTOPS, eRate[iBiexc]+hRate[iBiexc], 
              (eRate[iBiexc]+hRate[iBiexc])*AUTOPS, 1.0/((eRate[iBiexc]+hRate[iBiexc])*AUTOPS));      
            iBiexc++; // increments biexciton index energy range of biexciton was satisfied
          }
        }
      }
    }
  }
  fclose(pf);

  for (iDeltaE = 0; iDeltaE < numConsWindows; iDeltaE++) {
    writeSeparation(stdout);
    printf("Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calc_boltzmann_weighted_rates(biexcEnergy, &(eRate[iDeltaE*ist.numBiexcitons]), &(hRate[iDeltaE*ist.numBiexcitons]), par.temp, ist.numBiexcitons);
  }

  // Free dynamically allocated memory
  free(lista); free(listi);
  free(eRate); free(hRate); free(eCount); free(hCount); 
  free(biexcEnergy);

  return;
}

/*****************************************************************************/
// calculates and prints the boltzmann weighted AR, k_x-, k_x+ rates and 
// the corresponding lifetimes and boltzmann populations

void calc_boltzmann_weighted_rates(double *energies, double *eRate, double *hRate, double temp, int numStates) {
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
  printf("Boltzmann weighted k_x- = %.12f a.u. %.6f 1/ps\n", bwElecRate, bwElecRate*AUTOPS);
  printf("Boltzmann weighted k_x+ = %.12f a.u. %.6f 1/ps\n", bwHoleRate, bwHoleRate*AUTOPS);
  printf("Boltzmann weighted kAR  = %.12f a.u. %.6f 1/ps\n", bwAugerRate, bwAugerRate*AUTOPS);
  printf("Boltzmann weighted x- lifetime = %.6f a.u. %.12f ps\n", 1.0/(bwElecRate+EPS), 1.0/(bwElecRate*AUTOPS+EPS));
  printf("Boltzmann weighted x+ lifetime = %.6f a.u. %.12f ps\n", 1.0/(bwHoleRate+EPS), 1.0/(bwHoleRate*AUTOPS+EPS));
  printf("Boltzmann weighted Auger lifetime = %.6f a.u. %.12f ps\n", 1.0/(bwAugerRate+EPS), 1.0/(bwAugerRate*AUTOPS+EPS));

  return;
}

/*****************************************************************************/
// returns the partition function given an array of energies, a temperature,
// and the number of states (i.e., length of the energies array)

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
