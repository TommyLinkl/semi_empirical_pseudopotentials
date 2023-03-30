#include "fd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[],par_st *par,long_st *ist)
{
  FILE* pf = fopen("input.par" , "r");
  fscanf (pf,"%ld",&ist->nx);  /*** number of grid point in x ***/
  fscanf (pf,"%ld",&ist->ny);  /*** number of grid point in y ***/
  fscanf (pf,"%ld",&ist->nz);  /*** number of grid point in z ***/
  fscanf (pf,"%lg",&par->dx); /*** dx = minimum grid density ***/
  fscanf (pf,"%ld",&ist->ms);  /*** number of states per filter ***/
  fscanf (pf,"%ld",&ist->ns);  /*** number of filter cycles ***/
  fscanf (pf,"%ld",&ist->nc);    /*** length of newton interpolation ***/
  fscanf (pf,"%lg",&par->Elmin); /***minimum energy ***/
  fscanf (pf,"%lg",&par->Elmax); /*** maximum energy ***/
  fscanf (pf,"%ld",&ist->nthreads);
  fscanf (pf,"%ld",&ist->flagkb);
  fclose(pf);

  // Paramters that rarely, if ever, change - John
  ist->npot = 8192;
  par->Ekinmax = 20.0;

  // Get the number of atoms 
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);
  
  par->dy = par->dz = par->dx;
 
  return;
}

/*****************************************************************************/

void init_conf(double *rx,double *ry,double *rz,atm_st *atm,par_st *par,long_st *ist)
{
  long ntot, ntmp; FILE *pf;
  double xd, yd, zd;
  
  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  ist->nsemicond = read_conf(rx,ry,rz,atm,ntot,pf);
  fclose(pf);
  printf("nsemiconf = %ld\n",ist->nsemicond);
  
  xd = rint(0.5 * get_dot_ligand_size_z(rx, ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry, ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz, ntot) + 5.0);
  printf("xd = %g yd = %g zd = %g\n",xd,yd,zd);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  ntmp  = (long)((par->xmax - par->xmin) / par->dx);
  if (ntmp > ist->nx) ist->nx = ntmp;
  par->xmin = -((double)(ist->nx) * par->dx) / 2.0;
  par->xmax = ((double)(ist->nx) * par->dx) / 2.0;
  /*par->dx  = (par->xmax - par->xmin) / (double)(ist->nx);*/
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  ntmp  = (long)((par->ymax - par->ymin) / par->dy);
  if (ntmp > ist->ny) ist->ny = ntmp;
  /*par->dy  = (par->ymax - par->ymin) / (double)(ist->ny);*/
  par->ymin = -((double)(ist->ny) * par->dy) / 2.0;
  par->ymax = ((double)(ist->ny) * par->dy) / 2.0;
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  ntmp  = (long)((par->zmax - par->zmin) / par->dz) - 1;
  if (ntmp > ist->nz) ist->nz = ntmp;
  /*par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);*/
  par->zmin = -((double)(ist->nz) * par->dz) / 2.0;
  par->zmax = ((double)(ist->nz) * par->dz) / 2.0;
  par->dkz = TWOPI / ((double)ist->nz * par->dz);
  printf ("new xd = %g yd = %g zd = %g\n",par->xmax,par->ymax,par->zmax);

  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;

  ist->nspin = 2;
  ist->natomtype = 15;
  ist->mstot = ist->ms * ist->ns;
  ist->mstotngrid = ist->mstot * ist->ngrid;
  ist->nspinngrid = ist->nspin * ist->ngrid;

  printf ("nx = %ld  ny = %ld  nz = %ld npot = %ld\n",
	  ist->nx,ist->ny,ist->nz,ist->npot);
  printf ("ms = %ld  ns = %ld natom = %ld nc = %ld\n",
	  ist->ms,ist->ns,ist->natom,ist->nc);
  printf ("Elmin = %g  Elmax = %g  Ekinmax = %g\n",
	  par->Elmin,par->Elmax,par->Ekinmax);
  printf ("ngrid = %ld nspin = %ld\n",ist->ngrid,ist->nspin);
  printf ("threads = %ld\n",ist->nthreads);

  par->sigma = par->dx;
  par->sigma_1 = sqrt(0.5) / par->sigma;
  //par->Rnlcut2 = 0.49 * 6.0 * log(10.0) + 3.0 * par->sigma;
  par->Rnlcut2 = 0.49 * 6.0 * log(10.0);
  printf ("Rnlcut = %g\n",sqrt(par->Rnlcut2));
  fflush(stdout);

  return;
}

/*****************************************************************************/

void init_list(nlc_st *nlc,long *nl,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist)
{
  FILE *pf;
  long jatom, jx, jy, jz, jyz, jxyz;
  double dx, dy, dz, dxeps, dyeps, dzeps, dr_1, dr2;
  
  // Useful constants
  double tmp1 = 0.5 * sqrt(3.0 / PIE);
  double tmp2 = 0.5 * sqrt(3.0 / TWOPI);
  double tmp49 = 1.0 / 0.49;

  // Find all the grid points within par.Rnlcut of each atom and calculate
  // r, r2, Vso(r), etc. at the grid points and store the results in nlc
  for (jatom = 0; jatom < ist.nsemicond; jatom++){
    nl[jatom] = 0;
    for (jz = 0; jz < ist.nz; jz++){
      dz = vz[jz] - rz[jatom];
      dzeps = dz + EPSDX;
      for (jy = 0; jy < ist.ny; jy++){
      	jyz = ist.nx * (ist.ny * jz + jy);
      	dy = vy[jy] - ry[jatom];
      	dyeps = dy + EPSDX;
      	for (jx = 0; jx < ist.nx; jx++){
      	  jxyz = jyz + jx;
      	  dx = vx[jx] - rx[jatom];
      	  dxeps = dx + EPSDX;
      	  dr2 = dx * dx + dy * dy + dz * dz;
      	  if (dr2 < par.Rnlcut2){
      	    nlc[jatom*ist.nnlc + nl[jatom]].jxyz = jxyz;

      	    nlc[jatom*ist.nnlc + nl[jatom]].Vr = atm[jatom].Vso * exp(-tmp49 * dr2);
      	    nlc[jatom*ist.nnlc + nl[jatom]].r  = sqrt(dr2);

      	    dr_1 = 1.0 / sqrt(dx * dx + dy * dy + dzeps * dzeps);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y10.re = tmp1 * dzeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y10.im = 0.0;

      	    dr_1 = 1.0 / sqrt(dxeps * dxeps + dy * dy + dz * dz);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y11.re  = -tmp2 * dxeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1_1.re = tmp2 * dxeps * dr_1;

      	    dr_1 = 1.0 / sqrt(dx * dx + dyeps * dyeps + dz * dz);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y11.im  = -tmp2 * dyeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1_1.im = -tmp2 * dyeps * dr_1;

      	    if (dr2 > EPSDX) {
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2_1 = sqr(dr_1);
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2 = dr2;
      	    }
      	    else {
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2_1 = 0.0;
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2 = 1.0 / EPSDX;
      	    }
      	    nl[jatom]++;
      	  }
      	}
      }
    }
  }

  pf = fopen("list.dat" , "w");
  for (jatom = 0; jatom < ist.nsemicond; jatom++) {
    fprintf(pf, "%ld %ld\n", jatom, nl[jatom]);
  }
  fclose(pf);
  
  return;
}

/*****************************************************************************/

void init_kb(nlc_st *nlc,long *nl,double *Elkb,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist,double *dr,double *vr,double *potatom,long *npot)
{
  FILE *pf; char str[200], atype[3];
  long jatomtype, jtmp, jatom, jx, jy, jz, jgrid, index, natomtype, list[200];
  double dx, dy, dz, dr2, sum;
  double *ur1, *ur_2, rgrid, rgrid_1, *potPP;

  natomtype = get_number_of_atom_types(atm, ist, &list[0]);
  printf("number of atom types = %ld\n", natomtype);
  fflush(stdout);

  if ((ur1 = (double *) calloc(ist.npot*ist.natomtype, sizeof(double))) == NULL) nerror("ur1");
  if ((ur_2 = (double *) calloc(ist.npot*ist.natomtype, sizeof(double))) == NULL) nerror("ur_2");
  if ((potPP = (double *) calloc(ist.npot, sizeof(double)))==NULL) nerror("potPP");

  for (jatomtype = 0; jatomtype < natomtype; jatomtype++) {
    jatom = list[jatomtype+4]; // first 4 contain the index of the passivation ligands P1, P2, P3 and P4
    jtmp = jatom*ist.npot;
    assign_atom_type(atype, jatom);

    for (jgrid = 0; jgrid < ist.npot; jgrid++) {
      potPP[jgrid] = potatom[jtmp + jgrid];
    }

    if ((jatom==5) || (jatom==7)) {
      sprintf(str, "ur1%c.dat", atype[0]);
    }
    else if (jatom <= 12) {
      sprintf(str, "ur1%c%c.dat", atype[0], atype[1]);
    }
    else {
      sprintf(str, "ur1%c%c%c.dat", atype[0], atype[1], atype[2]); 
    }
    solve_for_ur(&vr[jtmp], &ur1[jtmp], potPP, dr[atm[jatom].natyp], npot[atm[jatom].natyp], 1.0, atm[jatom].initE, str);

    if ((jatom==5) || (jatom==7)) {
      sprintf(str, "ur_2%c.dat", atype[0]);
    }
    else if (jatom <= 12) {
      sprintf(str, "ur_2%c%c.dat", atype[0], atype[1]);
    }
    else {
      sprintf(str, "ur_2%c%c%c.dat", atype[0], atype[1], atype[2]);
    }
    solve_for_ur(&vr[jtmp], &ur_2[jtmp], potPP, dr[atm[jatom].natyp], npot[atm[jatom].natyp], -2.0, atm[jatom].initE, str);
    printf("%s\n", str);
    fflush(stdout);
  }
 
  for (jatom = 0; jatom < ist.nsemicond; jatom++){
    for (jgrid = 0, jz = 0; jz < ist.nz; jz++){
      dz = vz[jz] - rz[jatom];
      for (jy = 0; jy < ist.ny; jy++){
      	dy = vy[jy] - ry[jatom];
      	for (jx = 0; jx < ist.nx; jx++){
      	  dx = vx[jx] - rx[jatom];
      	  dr2 = dx * dx + dy * dy + dz * dz;
      	  if (dr2 < par.Rnlcut2) {
      	    rgrid = sqrt(dr2 + EPSDX);
      	    rgrid_1 = 1.0 / rgrid;
      	    index = jatom * ist.nnlc + jgrid;
      	    
      	    nlc[index].Vs0G1 = rgrid_1 * interpolate(rgrid, dr[atm[jatom].natyp], vr, ur1, ist.npot, npot[atm[jatom].natyp], atm[jatom].natyp);
      	    nlc[index].Vs0G_2 = rgrid_1 * interpolate(rgrid, dr[atm[jatom].natyp], vr, ur_2, ist.npot, npot[atm[jatom].natyp], atm[jatom].natyp);

      	    nlc[index].Vs0G1 *= nlc[index].Vr;
      	    nlc[index].Vs0G_2 *= nlc[index].Vr;
      	    
      	    jgrid++;
      	  }
      	}
      }
    }
  }

  for (jatom = 0; jatom < ist.nsemicond; jatom++) {
    for (sum = 0.0, jgrid = 0; jgrid < nl[jatom]; jgrid++) {
      index = jatom * ist.nnlc + jgrid;
      sum += sqr(nlc[index].r * nlc[index].Vs0G1) * nlc[index].Vr;
    }
    Elkb[jatom] = sum*par.dv;

    for (sum = 0.0, jgrid = 0; jgrid < nl[jatom]; jgrid++) {
      index = jatom * ist.nnlc + jgrid;
      sum += sqr(nlc[index].r * nlc[index].Vs0G_2) * nlc[index].Vr;
    }
    Elkb[jatom + ist.nsemicond] = sum*par.dv;
  }

  pf = fopen("Elkb.dat" , "w");
  for (jatom = 0; jatom < ist.nsemicond; jatom++) {
    fprintf(pf, "%ld %g %g\n", jatom, Elkb[jatom], Elkb[jatom+ist.nsemicond]);
  }
  fclose(pf);

  free(potPP); free(ur1); free(ur_2);

  return;
}

/*****************************************************************************/
// Calculates radial part of atomic solution of Dirac's equation by mapping
// radial Dirac equation into form that can be used as input to the
// generalized Numerov's method (in num_int.c) is used to solve the differential equation of the form
// u"(r) + p(r)*u'(r) + q(r)*u(r) + s(r) = 0
// for us:  s(r) = 0 and u(r) = R(r)/r, where R(r) is radial part of the wavefunction 

void solve_for_ur(double *r, double *ur, double *pot, double dr, long nr, double kappa, double initE, char *str) {
  FILE *pf;
  long i, jr, jnodes, jrMid; 
  double norm, sumIn, sumOut, E, dE, tmpE, potMin, scale, M, Vprime, c, c_2;
  double urInPrime, urOutPrime;	// used in correcting eigenvalue, E, guess 
  double *urIn, *urOut; 	// ur = urIn for r < rMid and ur = urOut fo r > rMid
  double *p, *q, *s;		// used to call gen_numerov_solver
  
  urIn = (double *) calloc(nr, sizeof(double));		// used for outward integration of u(r)
  urOut = (double *) calloc(nr, sizeof(double)); 	// used for inward integration of u(r)
  p = (double *) calloc(nr, sizeof(double));     	// p(r) = (1/(2Mc^2))*Vpp'
  q = (double *) calloc(nr, sizeof(double)); 		// q(r) = -[2M(Vpp-E)+k(k+1)/r^2-k/(2Mc^2r)*Vpp']
  s = (double *) calloc(nr, sizeof(double)); 		// s(r) = 0 for radial Dirac equation

  // Useful constants
  c = 137.036; 
  c_2 = 1.0/(c*c);
  potMin = min_element(pot, nr);
  E = initE;
  tmpE = initE;
  jrMid = 20;  				// r index where derivatives of urIn and urOut will be calculated

  for (i = 0; i < 1000; i++ ) {
    //printf("Entered the ur solver for kappa = % .8f on iteration %d with energy = % .12f\n", kappa, i, E);
    for (jnodes = 0; jnodes < 1000; jnodes++) {
      // Initial conditions
      if (! i && ! jnodes) {
        printf("Starting to calculate %s\n", str);
        urIn[0] = 0.0;	         		// initial conditions for the outward integration
        urIn[1] = 0.001;				// initial conditions for the outward integration
        //urOut[nr-1] = exp((E*r[nr-1])); 	// initial conditions for the inward integration
        //urOut[nr-2] = exp((E*r[nr-2]));      	// asymptotic behavior used
        urOut[nr-1] = 0.0001;      
        urOut[nr-2] = 0.00013;
       }
      // Calculate p(r), q(r), and s(r) to be used in the generalized Numerov solver
      for (jr = 0; jr < nr; jr++) {
        M = 1.0 - 0.5*c_2*(pot[jr]-E); 	// M = m - 1/(2c^2)*(Vpp-E)
        if (jr > 0 && jr < nr-1) Vprime = ((pot[jr+1]-pot[jr-1]) / (2.0*dr)); // first derivative of Vpp
        else if (jr == 0) Vprime = (pot[jr+1]-pot[jr]) / dr;
        else Vprime = (pot[jr]-pot[jr-1]) / dr;
        p[jr] = -0.5*c_2*Vprime/M;
        // TODO: EPSR0 and EPSR02 are very important since starting integration at r=0
        if (jr < 2) q[jr] = (2.0*M*(pot[jr]-E) + kappa*(kappa+1.0)/(r[jr]*r[jr]+EPSR02) - kappa*0.5*c_2*Vprime/(M*(r[jr]+EPSR0)));
        else q[jr] = (2.0*M*(pot[jr]-E) + kappa*(kappa+1.0)/(r[jr]*r[jr]) - kappa*0.5*c_2*Vprime/(M*(r[jr])));
        s[jr] = 0.0; // not really needed
      } 
      // Integrate from jr=0 to jr=jrMid (r < rmid)
      gen_numerov_solver(urIn, p, q, s, r, urIn[0], urIn[1], 0, jrMid+2);
      // make sure no nodes in urIn from r=0 to r=rmid - if node found, decrease E until no nodes found
      if (count_nodes(urIn, jrMid)) {      
        E -= 0.01;  // lower the energy to find a solution that does not have nodes inside rMid
        if (jnodes == 999) {
          tmpE += 10.0;
          printf("Ending energy = %.12f had %ld nodes\n", E, count_nodes(urIn, jrMid));
          E = tmpE;
          printf("Couldn't find state with zero nodes, new energy = %.12f\n", E);
        }
        continue;
      }
      else {
        break;  // exit out of the for loop over jnodes
      }
    }    

    // Integrate from jr=nr-1 to jr=jrMid
    gen_numerov_solver(urOut, p, q, s, r, urOut[nr-1], urOut[nr-2], nr-1, jrMid-3); 

    // Scale urOut such that urOut[jrMid]=urIn[jrMid] 
    scale = urIn[jrMid] / urOut[jrMid];
    for (jr = 0; jr < nr; jr++) urOut[jr] *= scale; 
 
    // Normalize ur
    // TODO: double check normalization and cite paper where found
    // normalization condition taken from: 
    sumIn = 0.0; sumOut = 0.0;
    for (jr = 0; jr < jrMid; jr++) sumIn += urIn[jr]*urIn[jr];
    sumOut = 0.0;
    for (jr = jrMid; jr < nr; jr++) sumOut += urOut[jr]*urOut[jr]; 
    norm = (sumIn + sumOut) * (4.0 * PIE * dr);
    norm = 1.0 / sqrt(norm+EPSDX);
    for (jr = 0; jr < nr; jr++) urIn[jr] *= norm;
    for (jr = 0; jr < nr; jr++) urOut[jr] *= norm; 
 
    // Calculate energy change, dE, based on first derivatives of urIn and urOut at rMid
    // TODO: double check the energy correction
    // procedure taken from: Comput. Phys. Commun. 184, 1777 (2013)  
    dE = energy_corrector_pt(urIn, urOut, dr, jrMid);
    // Break loop if energy eigenvalue has converged (i.e., ur converged)
    if (fabs(dE) < EPSE) {
      printf("Found an energy eigenvalue! E = % .12f\n", E); 
      break;
    }
    else if (i < 100) E += 0.1*dE;
    else E += dE;
  }

  // Store results in ur 
  for (jr = 0; jr < jrMid; jr++) ur[jr] = urIn[jr];
  for (jr = jrMid; jr < nr; jr++) ur[jr] = urOut[jr];

  // Print results to ur*.dat files
  pf = fopen(str , "w");
  for (jr = 0; jr < nr; jr++) {
    fprintf (pf,"% .8f % .8f % .8f % .8f % .8f % .8f % .8f\n", r[jr], ur[jr], urIn[jr], urOut[jr], pot[jr], p[jr], q[jr]);
  }
  fclose(pf);

  free(s); free(p); free(q);
  free(urOut); free(urIn);

  return;
}

/*****************************************************************************/
// returns how the energy should be updated to find the energy eigenvalue
// when solving the radial Dirac equation using the shooting outward and
// inward method - taken from Comput. Phys. Commun. 184, 1777 (2013)

double energy_corrector_pt(double *urIn, double *urOut, double dr, int jrMid) {
	double dE, urInPrime, urOutPrime;
	
	urInPrime = (urIn[jrMid+1]-urIn[jrMid-1])/(2.0*dr);
	urOutPrime = (urOut[jrMid+1]-urOut[jrMid-1])/(2.0*dr);

	dE = urOut[jrMid]*(urOutPrime-urInPrime);

	return dE;
}


/*****************************************************************************/

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,double *dr,double *vr,double *potatom,long *npot)
{
  FILE *pf;
  long jx, jy, jz, jyz, jxyz, jatom, ntot, jp, jtmp, nn, flags=0;
  double del, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double sum;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  
  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf ("dx = %g dy = %g dz = %g dv = %g dr = %g\n",
	  par->dx,par->dy,par->dz,par->dv,par->dr);

  /***initializing the ksqr vectors ***/
  for (ksqrx[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
    ksqrx[jx] = (ksqrx[ist->nx-jx] = 0.5 * sqr((double)(jx) * par->dkx) *
		ist->nx_1 * ist->ny_1 * ist->nz_1);

  for (ksqry[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
    ksqry[jy] = (ksqry[ist->ny-jy] = 0.5 * sqr((double)(jy) * par->dky) *
		ist->ny_1 * ist->nx_1 * ist->nz_1);

  for (ksqrz[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
    ksqrz[jz] = (ksqrz[ist->nz-jz] = 0.5 * sqr((double)(jz) * par->dkz) *
		ist->nz_1 * ist->nx_1 * ist->ny_1);

  par->Ekinmax *= (ist->ny_1 * ist->nx_1 * ist->nz_1);
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
      jxyz = jyz + jx;
      ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
      if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
    }
  }
  free(ksqrx); free(ksqry);  free(ksqrz);
  
  /*** read pseudopotentials ***/
  read_pot(vr,potatom,npot,dr,atm,ist->npot,ist->natomtype);
  
  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,jatom)
  for (jz = 0; jz < ist->nz; jz++) {
    for (jy = 0; jy < ist->ny; jy++) {
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++) {
      	jxyz = jyz + jx;
      	for (sum = 0.0, jatom = 0; jatom < ist->natom; jatom++){
      	  dx = vx[jx] - rx[jatom];
      	  dy = vy[jy] - ry[jatom];
      	  dz = vz[jz] - rz[jatom];
      	  del = sqrt(dx * dx + dy * dy + dz * dz);
      	  sum += interpolate(del,dr[atm[jatom].natyp],vr,potatom,ist->npot,npot[atm[jatom].natyp],atm[jatom].natyp);
      	}
      	potl[jxyz] = sum;
      }
    }
  }

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;
  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }

  printf("dV = %g vmin = %g Vmax = %g\n", par->Vmax-par->Vmin, par->Vmin, par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (par->dx*par->dx) + 0.5 * sqr(PIE) / (par->dy*par->dy) + 0.5 * sqr(PIE) / (par->dz*par->dz);
  printf("dT = %g\n", par->dE);
  
  /*** setting the energy grid El ***/
  del = (par->Elmax - par->Elmin) / (double)(ist->ms-1);
  for (jx = 0; jx < ist->ms; jx++) {
    eval[jx] = par->Elmin + (double)(jx) * del;
  }
  for (pf = fopen("Egrid.dat","w"),jx = 0; jx < ist->ms; jx++) {
    fprintf(pf, "%g\n", eval[jx]);
  }
  fclose(pf);
  
  for (nn = jatom = 0; jatom < ist->natom; jatom++)
    if (((atm[jatom].natyp < 8) || (atm[jatom].natyp > 11)) && (atm[jatom].natyp % 2)) nn++;
  printf("nn = %ld\n", nn);
  ist->homo = 4*nn-1;
  ist->lumo = ist->homo+1;
  printf("homo = %ld lumo = %ld\n", ist->homo, ist->lumo);

  
  /*** for the nonlocal potential ***/
  for (ist->nnlc = 0, jatom = 0; jatom < ist->nsemicond; jatom++) {
    for (jtmp =0, jz = 0; jz < ist->nz; jz++) {
      for (jy = 0; jy < ist->ny; jy++) {
      	for (jx = 0; jx < ist->nx; jx++) {
      	  dx = vx[jx] - rx[jatom];
      	  dy = vy[jy] - ry[jatom];
      	  dz = vz[jz] - rz[jatom];
      	  if (dx*dx+dy*dy+dz*dz < par->Rnlcut2) {
            jtmp++; 
          }
      	}
      }
    }
    if (jtmp > ist->nnlc) ist->nnlc = jtmp;
  }
  printf("nnlc = %d\n", ist->nnlc);
  fflush(stdout);

  return;
}
/*****************************************************************************/

void init_psi(zomplex *psi,long_st ist,par_st par,long *idum)
{
  long jx, jy, jz, jzy, jxyz;
  long tidum = (*idum);

  for (jz = 0; jz < ist.nz; jz++) for (jy = 0; jy < ist.ny; jy++) {
    for (jzy = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
      jxyz = jzy + jx;
      psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
      psi[jxyz].im = 0.0;
    }
  }

  normalize(psi,par.dv,ist.ngrid,ist.nthreads);
  (*idum) = tidum;

  return;
}

/*****************************************************************************/
