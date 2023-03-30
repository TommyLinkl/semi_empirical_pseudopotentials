/*****************************************************************************/
// File contains 

#include "fd.h"

/*****************************************************************************/

void hamiltonian(zomplex *phi, zomplex *psi, double *potl, double *ksqr, long_st ist, par_st par, nlc_st *nlc,
                long *nl, double *Elkb, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi) {
  long i, itmp;
  zomplex *psinl;

  // Dynamically allocate memory 
  if ((psinl = (zomplex *) calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psinl");
  
  // Copy phi into psi
  memcpy(&psi[0], &phi[0], ist.nspinngrid*sizeof(psi[0]));

  // Calculate the action of the kinetic energy part of the Hamiltonian on psi: |phi> = T|psi>
  kinetic(&phi[0*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist);
  kinetic(&phi[1*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist);

  // Calculate the action of the non-local spin-orbit part of the Hamiltonian on psi: |psinl> = Vso|psi>
  if (ist.flagkb == 1) {
    nonlocal_potential_kb(psinl, psi, ist, par, nlc, nl, Elkb);
  }
  else {
    nonlocal_potential(psinl, psi, ist, par, nlc, nl);
  } 
    
  // Calculate the action of the local potential energy part of the Hamiltonian on psi
  // and add them all together: |phi> = T|psi> + Vlocal|psi> + Vso|psi>
  for (i = 0; i < ist.ngrid; i++) {
    phi[i].re += (potl[i] * psi[i].re + psinl[i].re);
    phi[i].im += (potl[i] * psi[i].im + psinl[i].im);

    itmp = i + ist.ngrid;
    phi[itmp].re += (potl[i] * psi[itmp].re + psinl[itmp].re);
    phi[itmp].im += (potl[i] * psi[itmp].im + psinl[itmp].im);
  }
  
  // Free dynamically allocated memory
  free(psinl);

  return;
}

/*****************************************************************************/
// Calculates T|psi> via FFT and stores result in psi 

void kinetic(zomplex *psi, double *ksqr, fftw_plan_loc planfw, fftw_plan_loc planbw,
              fftw_complex *fftwpsi,long_st ist)
{
  long i;

  // Copy psi to fftwpsi
  memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));
  
  // FT from r-space to k-space
  fftw_execute(planfw);
  
  // Kinetic energy is diagonal in k-space, just multiply fftwpsi by k^2
  for (i = 0; i < ist.ngrid; i++) {
    fftwpsi[i][0] *= ksqr[i];
    fftwpsi[i][1] *= ksqr[i];
  }

  // Inverse FT back to r-space
  fftw_execute(planbw);
  
  // Copy fftwpsi to psi to store T|psi> into |psi>
  memcpy(&psi[0], &fftwpsi[0], ist.ngrid*sizeof(psi[0]));

  return;
}

/*****************************************************************************/
// Calculates the action of the spin-orbit, nonlocal potential on the complex
// spinor and stores the result in psinl: |psinl> = Vso|psi>

void nonlocal_potential(zomplex *psinl, zomplex *psi, long_st ist, par_st par, nlc_st *nlc, long *nl)
{
  long jatom, j1, j2; 
  long index1, index2, index3, index1_psi, index2_psi, index3_psi;
  double tmpexp, r_rprime; 
  zomplex tmp10, tmp11, tmp1_1, sum11up, sum10up, sum1_1up, sum11dn, sum10dn, sum1_1dn;

  // Set useful constants
  double tmp = (0.5 / sqrt(TWOPI * sqr(par.sigma)) * par.dv); // 0.5 comes from Vso/2, dv for the intergration over r' 
  double sqrttwo = sqrt(2.0);

  // Zero psinl 
  for (j1 = 0; j1 < ist.ngrid; j1++) {
    psinl[j1].re = psinl[j1].im = 0.0;
  }
  
  // Calculate the action of Vso for each atom as they act independently 
  for (jatom = 0; jatom < ist.nsemicond; jatom++) {
    // Loop over all grid points (nl[jatom]) in the neighborhood of jatom that were previously 
    // stored in nlc[jatom*ist.nnlc] and calculate psinl at each of these grid points
    for (j1 = 0; j1 < nl[jatom]; j1++) {
      index1 = jatom * ist.nnlc + j1;

      // The index of the (x,y,z,up) grid point that we are calculating psinl at (i.e., <r|psinl>)
      index1_psi = nlc[index1].jxyz;  

      // The part that does not depend on r' and is, thus, outside of the integral of j2          
      tmp10.re = tmp * nlc[index1].y10.re;
      tmp10.im = tmp * nlc[index1].y10.im;

      tmp11.re = tmp * nlc[index1].y11.re;
      tmp11.im = tmp * nlc[index1].y11.im;

      tmp1_1.re = tmp * nlc[index1].y1_1.re;
      tmp1_1.im = tmp * nlc[index1].y1_1.im;

      // Zero the parts that will store the result of the integration
      sum11up.re = sum11up.im = sum10up.re = sum10up.im = sum1_1up.re = sum1_1up.im = 0.0;
      sum11dn.re = sum11dn.im = sum10dn.re = sum10dn.im = sum1_1dn.re = sum1_1dn.im = 0.0;

      // Here is the non-local part of the spin-orbit pseudopotentials: 
      // we have to integrate over all the grid points near jatom (indexed by j2) 
      // in order to obtain <r|psinl> = <r|Vso|psi> where r is indexed by j1 
      for (j2 = 0; j2 < nl[jatom]; j2++) {
      	index2 = jatom * ist.nnlc + j2;

        // r-r' to be used in Gaussain representation of delta(r-r')
      	r_rprime = nlc[index1].r - nlc[index2].r; 

        // The 0.5*(nlc[index1].Vr + nlc[index2].Vr) is just the average of Vso(r) as index1 and index2 
        // should be the same distance from the atom and, thus, equal (as required by the delta(r-r'))
        // the Gaussian is used to represent delta(r-r') and sigma_1 = 1/(sqrt(2.0)*sigma)
        // we then divide by (r')^2 where we have agained used the averaging trick
        tmpexp = 0.5 * (nlc[index1].Vr + nlc[index2].Vr) *
                  exp(-sqr(r_rprime * par.sigma_1)) / (0.5 * (nlc[index1].r2 + nlc[index2].r2));

        // index2_psi is used to obtain the spin-up component of psi at r' ( <r'up|psi> )
      	index2_psi = nlc[index2].jxyz;  

      	sum11up.re +=  tmpexp * (nlc[index2].y11.re * psi[index2_psi].re + nlc[index2].y11.im * psi[index2_psi].im);
      	sum11up.im +=  tmpexp * (nlc[index2].y11.re * psi[index2_psi].im - nlc[index2].y11.im * psi[index2_psi].re);
      	
      	sum10up.re +=  tmpexp * (nlc[index2].y10.re * psi[index2_psi].re + nlc[index2].y10.im * psi[index2_psi].im);
      	sum10up.im +=  tmpexp * (nlc[index2].y10.re * psi[index2_psi].im - nlc[index2].y10.im * psi[index2_psi].re);
      	
      	sum1_1up.re +=  tmpexp * (nlc[index2].y1_1.re * psi[index2_psi].re + nlc[index2].y1_1.im * psi[index2_psi].im);
      	sum1_1up.im +=  tmpexp * (nlc[index2].y1_1.re * psi[index2_psi].im - nlc[index2].y1_1.im * psi[index2_psi].re);
      	
        // index3_psi is used to obtain the spin-down component of psi at r' ( <r'down|psi> )
      	index3_psi = ist.ngrid + index2_psi;

      	sum11dn.re +=  tmpexp * (nlc[index2].y11.re * psi[index3_psi].re + nlc[index2].y11.im * psi[index3_psi].im);
      	sum11dn.im +=  tmpexp * (nlc[index2].y11.re * psi[index3_psi].im - nlc[index2].y11.im * psi[index3_psi].re);
      	
      	sum10dn.re +=  tmpexp * (nlc[index2].y10.re * psi[index3_psi].re + nlc[index2].y10.im * psi[index3_psi].im);
      	sum10dn.im +=  tmpexp * (nlc[index2].y10.re * psi[index3_psi].im - nlc[index2].y10.im * psi[index3_psi].re);
      	
      	sum1_1dn.re +=  tmpexp * (nlc[index2].y1_1.re * psi[index3_psi].re + nlc[index2].y1_1.im * psi[index3_psi].im);
      	sum1_1dn.im +=  tmpexp * (nlc[index2].y1_1.re * psi[index3_psi].im - nlc[index2].y1_1.im * psi[index3_psi].re);
      }
      
      // for spin-up component of the spinor
      psinl[index1_psi].re += (tmp11.re * sum11up.re - tmp11.im * sum11up.im) -
                                (tmp1_1.re * sum1_1up.re - tmp1_1.im * sum1_1up.im) +
                                sqrttwo * ((tmp10.re * sum11dn.re - tmp10.im * sum11dn.im) +
                                (tmp1_1.re * sum10dn.re - tmp1_1.im * sum10dn.im));
      psinl[index1_psi].im += (tmp11.im * sum11up.re + tmp11.re * sum11up.im) -
                                (tmp1_1.im * sum1_1up.re + tmp1_1.re * sum1_1up.im) +
                                sqrttwo * ((tmp10.im * sum11dn.re + tmp10.re * sum11dn.im) +
                                (tmp1_1.im * sum10dn.re + tmp1_1.re * sum10dn.im));
      
      // for spin-down component of the spinor
      index3_psi = ist.ngrid + index1_psi;
      psinl[index3_psi].re += (-tmp11.re * sum11dn.re + tmp11.im * sum11dn.im) -
                                (tmp1_1.re * sum1_1dn.re - tmp1_1.im * sum1_1dn.im) +
                                sqrttwo * ((tmp10.re * sum1_1up.re - tmp10.im * sum1_1up.im) +
                                (tmp11.re * sum10up.re - tmp11.im * sum10up.im));
      psinl[index3_psi].im += (-tmp11.im * sum11dn.re - tmp11.re * sum11dn.im) -
                                (tmp1_1.im * sum1_1dn.re + tmp1_1.re * sum1_1dn.im) +
                                sqrttwo * ((tmp10.im * sum1_1up.re + tmp10.re * sum1_1up.im) +
                                (tmp11.im * sum10up.re + tmp11.re * sum10up.im));
    }
  }

  return;
}

/*****************************************************************************/
// Calculates the action of the nonlocal, spin-orbit potential on the complex
// spinor ket using the Kleinman-Bylander technique and stores the result in psinl
// |psinl> = Vso|psi> 

void nonlocal_potential_kb(zomplex *psinl, zomplex *psi, long_st ist, par_st par, 
                            nlc_st *nlc, long *nl, double *Elkb) {
  ing_st tmping; 
  zomplex tmpsi;
  long jatom, jgrid, index, index_psi;
  double tmpVs0, tmpVs0G1, tmpVs0G_2, tmp;
  double ElG1, ElG_2;
  
  // Set useful constants
  double oneoversix = 1.0 / 6.0;
  double oneoverthree = 1.0 / 3.0;
  double twooverthree = 2.0 / 3.0;
  double sqrttwooversix = sqrt(2.0) / 6.0;
  double sqrttwooverthree = sqrt(2.0) / 3.0;

  // Zero psinl
  for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
    psinl[jgrid].re = psinl[jgrid].im = 0.0;
  }
  
  // Calculate the action of Vso for each atom as they act independently
  for (jatom = 0; jatom < ist.nsemicond; jatom++) {

    // Zero the parts that will store the result of the integration
    tmping.G1y10up.re = tmping.G1y10up.im = 0.0;
    tmping.G1y11up.re = tmping.G1y11up.im = 0.0;
    tmping.G1y1_1up.re = tmping.G1y1_1up.im = 0.0;
    tmping.G1y10dn.re = tmping.G1y10dn.im = 0.0;
    tmping.G1y11dn.re = tmping.G1y11dn.im = 0.0;
    tmping.G1y1_1dn.re = tmping.G1y1_1dn.im = 0.0;

    tmping.G_2y10up.re = tmping.G_2y10up.im = 0.0;
    tmping.G_2y11up.re = tmping.G_2y11up.im = 0.0;
    tmping.G_2y1_1up.re = tmping.G_2y1_1up.im = 0.0;
    tmping.G_2y10dn.re = tmping.G_2y10dn.im = 0.0;
    tmping.G_2y11dn.re = tmping.G_2y11dn.im = 0.0;
    tmping.G_2y1_1dn.re = tmping.G_2y1_1dn.im = 0.0;

    for (jgrid = 0; jgrid < nl[jatom]; jgrid++) {
      index = jatom * ist.nnlc + jgrid;

      // The index of the (x,y,z,up) grid point that we are calculating psinl at (i.e., <r|psinl>)
      index_psi = nlc[index].jxyz;

      // Value of the spin-up component of psi at r, <r|psi>
      tmpsi = psi[index_psi];

      tmpVs0 = nlc[index].Vs0G1;
      tmping.G1y10up.re += tmpVs0 * (nlc[index].y10.re * tmpsi.re + nlc[index].y10.im * tmpsi.im);
      tmping.G1y10up.im += tmpVs0 * (nlc[index].y10.re * tmpsi.im - nlc[index].y10.im * tmpsi.re);
      tmping.G1y11up.re += tmpVs0 * (nlc[index].y11.re * tmpsi.re + nlc[index].y11.im * tmpsi.im);
      tmping.G1y11up.im += tmpVs0 * (nlc[index].y11.re * tmpsi.im - nlc[index].y11.im * tmpsi.re);
      tmping.G1y1_1up.re += tmpVs0 * (nlc[index].y1_1.re * tmpsi.re + nlc[index].y1_1.im * tmpsi.im);
      tmping.G1y1_1up.im += tmpVs0 * (nlc[index].y1_1.re * tmpsi.im - nlc[index].y1_1.im * tmpsi.re);

      tmpVs0 = nlc[index].Vs0G_2;
      tmping.G_2y10up.re += tmpVs0 * (nlc[index].y10.re * tmpsi.re + nlc[index].y10.im * tmpsi.im);
      tmping.G_2y10up.im += tmpVs0 * (nlc[index].y10.re * tmpsi.im - nlc[index].y10.im * tmpsi.re);
      tmping.G_2y11up.re += tmpVs0 * (nlc[index].y11.re * tmpsi.re + nlc[index].y11.im * tmpsi.im);
      tmping.G_2y11up.im += tmpVs0 * (nlc[index].y11.re * tmpsi.im - nlc[index].y11.im * tmpsi.re);
      tmping.G_2y1_1up.re += tmpVs0 * (nlc[index].y1_1.re * tmpsi.re + nlc[index].y1_1.im * tmpsi.im);
      tmping.G_2y1_1up.im += tmpVs0 * (nlc[index].y1_1.re * tmpsi.im - nlc[index].y1_1.im * tmpsi.re);
      
      // The index of the (x,y,z,down) grid point that we are calculating psinl at (i.e., <r|psinl>)
      index_psi = ist.ngrid + index_psi;

      // Value of the spin-down component of psi at r, <r|psi>
      tmpsi = psi[index_psi];

      tmpVs0 = nlc[index].Vs0G1;
      tmping.G1y10dn.re += tmpVs0 * (nlc[index].y10.re * tmpsi.re + nlc[index].y10.im * tmpsi.im);
      tmping.G1y10dn.im += tmpVs0 * (nlc[index].y10.re * tmpsi.im - nlc[index].y10.im * tmpsi.re);
      tmping.G1y11dn.re += tmpVs0 * (nlc[index].y11.re * tmpsi.re + nlc[index].y11.im * tmpsi.im);
      tmping.G1y11dn.im += tmpVs0 * (nlc[index].y11.re * tmpsi.im - nlc[index].y11.im * tmpsi.re);
      tmping.G1y1_1dn.re += tmpVs0 * (nlc[index].y1_1.re * tmpsi.re + nlc[index].y1_1.im * tmpsi.im);
      tmping.G1y1_1dn.im += tmpVs0 * (nlc[index].y1_1.re * tmpsi.im - nlc[index].y1_1.im * tmpsi.re);

      tmpVs0 = nlc[index].Vs0G_2;
      tmping.G_2y10dn.re += tmpVs0 * (nlc[index].y10.re * tmpsi.re + nlc[index].y10.im * tmpsi.im);
      tmping.G_2y10dn.im += tmpVs0 * (nlc[index].y10.re * tmpsi.im - nlc[index].y10.im * tmpsi.re);
      tmping.G_2y11dn.re += tmpVs0 * (nlc[index].y11.re * tmpsi.re + nlc[index].y11.im * tmpsi.im);
      tmping.G_2y11dn.im += tmpVs0 * (nlc[index].y11.re * tmpsi.im - nlc[index].y11.im * tmpsi.re);
      tmping.G_2y1_1dn.re += tmpVs0 * (nlc[index].y1_1.re * tmpsi.re + nlc[index].y1_1.im * tmpsi.im);
      tmping.G_2y1_1dn.im += tmpVs0 * (nlc[index].y1_1.re * tmpsi.im - nlc[index].y1_1.im * tmpsi.re);
    }

    // TODO: check this part, J.P. May 2nd 2019
    // 1) Is par.dv needed? Elkb[..] already got multiplied by par.dv in init_kb
    // 2) Should it be 1/(Elkb[..])? 
    ElG1 = Elkb[jatom] * par.dv;
    ElG_2 = Elkb[jatom+ist.nsemicond] * par.dv;
        
    tmping.G1y10up.re *= ElG1;
    tmping.G1y10up.im *= ElG1; 
    tmping.G1y11up.re *= ElG1;
    tmping.G1y11up.im *= ElG1; 
    tmping.G1y1_1up.re *= ElG1;
    tmping.G1y1_1up.im *= ElG1; 

    tmping.G1y10dn.re *= ElG1;
    tmping.G1y10dn.im *= ElG1; 
    tmping.G1y11dn.re *= ElG1;
    tmping.G1y11dn.im *= ElG1; 
    tmping.G1y1_1dn.re *= ElG1;
    tmping.G1y1_1dn.im *= ElG1; 

    tmping.G_2y10up.re *= ElG_2;
    tmping.G_2y10up.im *= ElG_2; 
    tmping.G_2y11up.re *= ElG_2;
    tmping.G_2y11up.im *= ElG_2; 
    tmping.G_2y1_1up.re *= ElG_2;
    tmping.G_2y1_1up.im *= ElG_2; 

    tmping.G_2y10dn.re *= ElG_2;
    tmping.G_2y10dn.im *= ElG_2; 
    tmping.G_2y11dn.re *= ElG_2;
    tmping.G_2y11dn.im *= ElG_2; 
    tmping.G_2y1_1dn.re *= ElG_2;
    tmping.G_2y1_1dn.im *= ElG_2; 

    for (jgrid = 0; jgrid < nl[jatom]; jgrid++) {
      index = jatom * ist.nnlc + jgrid; 

      // The index of the (x,y,z,up) grid point that we are calculating psinl at (i.e., <r|psinl>)
      index_psi = nlc[index].jxyz; 

      tmpVs0G1 = nlc[index].Vs0G1;  
      tmpVs0G_2 = nlc[index].Vs0G_2;  

      tmp = oneoversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y1_1up.re * nlc[index].y1_1.re - tmping.G_2y1_1up.im * nlc[index].y1_1.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y1_1up.re * nlc[index].y1_1.im + tmping.G_2y1_1up.im * nlc[index].y1_1.re);

      tmp = oneoverthree * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y10up.re * nlc[index].y10.re - tmping.G_2y10up.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y10up.re * nlc[index].y10.im + tmping.G_2y10up.im * nlc[index].y10.re);
      
      tmp = 0.5 * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y11up.re * nlc[index].y11.re - tmping.G_2y11up.im * nlc[index].y11.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y11up.re * nlc[index].y11.im + tmping.G_2y11up.im * nlc[index].y11.re);

      tmp = sqrttwooversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y10dn.re * nlc[index].y1_1.re - tmping.G_2y10dn.im * nlc[index].y1_1.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y10dn.re * nlc[index].y1_1.im + tmping.G_2y10dn.im * nlc[index].y1_1.re);

      tmp = sqrttwooversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y11dn.re * nlc[index].y10.re - tmping.G_2y11dn.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y11dn.re * nlc[index].y10.im + tmping.G_2y11dn.im * nlc[index].y10.re);
      
      tmp = -twooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y1_1up.re * nlc[index].y1_1.re - tmping.G1y1_1up.im * nlc[index].y1_1.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y1_1up.re * nlc[index].y1_1.im + tmping.G1y1_1up.im * nlc[index].y1_1.re);

      tmp = -oneoverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y10up.re * nlc[index].y10.re - tmping.G1y10up.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y10up.re * nlc[index].y10.im + tmping.G1y10up.im * nlc[index].y10.re);

      tmp = sqrttwooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y10dn.re * nlc[index].y1_1.re - tmping.G1y10dn.im * nlc[index].y1_1.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y10dn.re * nlc[index].y1_1.im + tmping.G1y10dn.im * nlc[index].y1_1.re);

      tmp = sqrttwooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y11dn.re * nlc[index].y10.re - tmping.G1y11dn.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y11dn.re * nlc[index].y10.im + tmping.G1y11dn.im * nlc[index].y10.re);

      // The index of the (x,y,z,down) grid point that we are calculating psinl at (i.e., <r|psinl>)
      index_psi = ist.ngrid + index_psi; 

      tmp = oneoversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y11dn.re * nlc[index].y11.re - tmping.G_2y11dn.im * nlc[index].y11.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y11dn.re * nlc[index].y11.im + tmping.G_2y11dn.im * nlc[index].y11.re);

      tmp = oneoverthree * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y10dn.re * nlc[index].y10.re - tmping.G_2y10dn.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y10dn.re * nlc[index].y10.im + tmping.G_2y10dn.im * nlc[index].y10.re);

      tmp = -0.5 * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y1_1dn.re * nlc[index].y1_1.re - tmping.G_2y1_1dn.im * nlc[index].y1_1.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y1_1dn.re * nlc[index].y1_1.im + tmping.G_2y1_1dn.im * nlc[index].y1_1.re);

      tmp = sqrttwooversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y1_1up.re * nlc[index].y10.re - tmping.G_2y1_1up.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y1_1up.re * nlc[index].y10.im + tmping.G_2y1_1up.im * nlc[index].y10.re);

      tmp = sqrttwooversix * tmpVs0G_2;
      psinl[index_psi].re +=  tmp * (tmping.G_2y10up.re * nlc[index].y11.re - tmping.G_2y10up.im * nlc[index].y11.im);
      psinl[index_psi].im +=  tmp * (tmping.G_2y10up.re * nlc[index].y11.im + tmping.G_2y10up.im * nlc[index].y11.re);

      tmp = sqrttwooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y1_1up.re * nlc[index].y10.re - tmping.G1y1_1up.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y1_1up.re * nlc[index].y10.im + tmping.G1y1_1up.im * nlc[index].y10.re);

      tmp = sqrttwooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y10up.re * nlc[index].y10.re - tmping.G1y10up.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y10up.re * nlc[index].y10.im + tmping.G1y10up.im * nlc[index].y10.re);

      tmp = -oneoverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y10dn.re * nlc[index].y10.re - tmping.G1y10dn.im * nlc[index].y10.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y10dn.re * nlc[index].y10.im + tmping.G1y10dn.im * nlc[index].y10.re);

      tmp = -twooverthree * tmpVs0G1;
      psinl[index_psi].re +=  tmp * (tmping.G1y11dn.re * nlc[index].y11.re - tmping.G1y11dn.im * nlc[index].y11.im);
      psinl[index_psi].im +=  tmp * (tmping.G1y11dn.re * nlc[index].y11.im + tmping.G1y11dn.im * nlc[index].y11.re);
    }
  }
  
  return;
}

/*****************************************************************************/