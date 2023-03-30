/****************************************************************************/

#include "nc.h"
#include "unistd.h"

/****************************************************************************/

void readInput(param *params) {
	FILE *pf;
	int j, i = 0;
	char field[1000], tmp[1000];

	// Set defaults
	strcpy(params->confFile, "conf.par");
	strcpy(params->confFileUnits, "atomic");
	strcpy(params->centerAtomSymbol, "Cd"); 
	strcpy(params->ncType, "QD");
	strcpy(params->mainNC.ncCrystalStructure, "wurtzite");
	params->nMaxBonds = 4;
	params->nMinBonds = params->nMaxBonds-2;
	params->remDanglingAtoms = 0;
	params->nLayersToCut = 0;
	params->passivate = 0;
	params->makeChiralNC = 0;
	params->buildNewNC = 0;
	params->buildNewBulk = 1; // used if buildNewNC equals 1 (true) in input.par
	params->nStackingFaults = 0;
	params->nCores = 0;
	params->nAttachedNCs = 0;
	strcpy(params->mainNC.ncAtomSymbol1, "Cd"); strcpy(params->mainNC.ncAtomSymbol1, "Se"); 

	// Open and read input.par if it exists - use defaults otherwise 
	if( access( "input.par", F_OK) != -1 ) {
		pf = fopen("input.par", "r");
		while (fscanf(pf, "%s", field) != EOF && i < 13) {
			if (! strcmp(field, "confFile")) {
				fscanf(pf, "%s %s", tmp, &(params->confFile));
				params->buildNewBulk = 0;
			}
			else if (! strcmp(field, "confFileUnits")) fscanf(pf, "%s %s", tmp, &(params->confFileUnits));
			else if (! strcmp(field, "nMaxBonds")) fscanf(pf, "%s %d", tmp, &(params->nMaxBonds));
			else if (! strcmp(field, "remDanglingAtoms")) fscanf(pf, "%s %d", tmp, &(params->remDanglingAtoms));
			else if (! strcmp(field, "passivate")) fscanf(pf, "%s %d", tmp, &(params->passivate));
			else if (! strcmp(field, "centerAtomSymbol")) fscanf(pf, "%s %s", tmp, &(params->centerAtomSymbol));
			else if (! strcmp(field, "ncType"))  fscanf(pf, "%s %s", tmp, &(params->ncType));
			else if (! strcmp(field, "nLayersToCut")) {
				fscanf(pf, "%s %d", tmp, &(params->nLayersToCut));
				params->remDanglingAtoms = 1; // always used when cutting layers off NC
			} 
			else if (! strcmp(field, "buildNewNC")) {
				fscanf(pf, "%s %d", tmp, &(params->buildNewNC));
				if (params->buildNewNC) {
					params->mainNC.index = 0;
					fillNanocrystalStructure(&(params->mainNC), pf);
					strcpy(params->ncType, params->mainNC.ncType);
				}
			}
			else if (! strcmp(field, "nCores")) {
				fscanf(pf, "%s %d", tmp, &(params->nCores));
				for (j = 0; j < params->nCores; j++) {
					params->coreNC[j].index = j;
					fillNanocrystalStructure(&(params->coreNC[j]), pf);
				}
			}
			else if  (! strcmp(field, "nAttachedNCs")) {
				fscanf(pf, "%s %d", tmp, &(params->nAttachedNCs));
				for (j = 0; j < params->nAttachedNCs; j++) {
					params->attachedNC[j].index = j;
					fillNanocrystalStructure(&(params->attachedNC[j]), pf);
				}
			}
			else if (! strcmp(field, "nStackingFaults")) {
				fscanf(pf, "%s %d", tmp, &(params->nStackingFaults));
				for (j = 0; j < params->nStackingFaults; j++) {
					fscanf(pf, "%s %s %d", field, tmp, &(params->stackFault[j].type));
					fscanf(pf, "%s %s %d", field, tmp, &(params->stackFault[j].position));					
				}	
			}
			else if (! strcmp(field, "makeChiralNC")) {
				fscanf(pf, "%s %d", tmp, &(params->makeChiralNC));
				fscanf(pf, "%s %s %d", field, tmp, &(params->chiralNanocrystalParams.type));
				fscanf(pf, "%s %s %s", field, tmp, &(params->chiralNanocrystalParams.chiralAtomSymbol));					
			} 
			else {
				printf("Invalid input field and/ or format - equal sign required after each field\n");
				printf("Only allowed fields are (case-sensitive):\n\n");
				printf("confFile = fileName (nAtoms first line then atomSymbol x y z\n");
				printf("confFileUnits = atomic or angstroms\n");
				printf("centerAtomSymbol = Se (Cd is default - use XX for no center)\n");
				printf("nMaxBonds = 4 (CdSe/CdS/ZnS/..)\n");
				printf("buildNewNC = 1 (0 -> false is default - must be followed immediately by ncType and size)\n");
				printf("ncType = NPL (QD is default (QD, NR, NPL or attachedNCs) - must be below buildNewNC and nCores\n");
				printf("nCores = 1 (0 -> false is default (0,1,2,...) - bulk crystal built with lattice constant of the mainNC)\n\n");
				printf("nAttachedNCs = 2 (0 -> false is default (1,3) - bulk crystal built with lattice constant of the mainNC)\n\n");
				printf("nStackingFaults = 1, 2,... (0 -> false is default)\n");
				printf("makeChiralNC = 1.(0 -> false is default - followed by chiralType = int and chiralAtomSymbol = Cd (Se..) lines)\n");
				printf("remDanglingAtoms = 1 (0 -> false is default)\n");
				printf("passivate = 1 (0 -> false is default)\n");
				printf("nLayersToCut = 0, 1, 2,... (0 is default)\n\n"); fflush(stdout);
				exit(EXIT_FAILURE);
			}
			i++;
		}
		fclose(pf);
	}
	else {
		printf("No input.par file detected in cwd - using defaults!\n");
	}

	// Print input used to screen
	printf("INPUT PARAMETERS:\n\n");
	if (params->buildNewNC) {
		printf("buildNewNC = %d\n", params->buildNewNC);
		printf("centerAtomSymbol = %s\n", params->centerAtomSymbol);
		printf("nCores = %d\n", params->nCores);
		printf("nAttachedNCs = %d\n", params->nAttachedNCs);
		printf("ncType = %s\n", params->mainNC.ncType);
		printf("atomSymbols = %s %s\n", params->mainNC.ncAtomSymbol1, params->mainNC.ncAtomSymbol2);
		printf("ncCrystalStructure = %s\n", params->mainNC.ncCrystalStructure);
		if (! strcmp(params->mainNC.ncType, "QD")) {
			printf("qdDiameter = %.2f A\n", params->mainNC.ncSize.mag);
		}
		else if (! strcmp(params->mainNC.ncType, "NR")) {
			printf("nrDiameter = %.2f A\n", params->mainNC.ncSize.mag);
			printf("nrLength = %.2f A\n", params->mainNC.ncSize.z);
		} 
		else if (! strcmp(params->mainNC.ncType, "NPL")) {
			printf("nplLength = %.2f A\n", params->mainNC.ncSize.x);
			printf("nplWidth = %.2f A\n", params->mainNC.ncSize.y);
			printf("nplThickness = %.2f A\n", params->mainNC.ncSize.z);
		}
		for (i = 0; i < params->nCores; i++) {
			printf("Core %d is a %s-%s %s\n", i, params->coreNC[i].ncAtomSymbol1, 
					params->coreNC[i].ncAtomSymbol2, params->coreNC[i].ncType);
			printf("The origin of this core is = %.2f % .2f % .2f\n", params->coreNC[i].ncCenter.x,
						params->coreNC[i].ncCenter.y, params->coreNC[i].ncCenter.z);
			printf("The size of this core is = ");
			if (! strcmp(params->coreNC[i].ncType, "QD")) {
				printf("D = %.2f A\n", params->coreNC[i].ncSize.mag);
			}
			else if (! strcmp(params->coreNC[i].ncType, "NR")) {
				printf("D = %.2f A ", params->coreNC[i].ncSize.mag);
				printf("L = %.2f A\n", params->coreNC[i].ncSize.z);
			}
			else if (! strcmp(params->coreNC[i].ncType, "NPL")) {
				printf("L = %.2f A ", params->coreNC[i].ncSize.x);
				printf("W = %.2f A ", params->coreNC[i].ncSize.y);
				printf("H = %.2f A\n", params->coreNC[i].ncSize.z);
			}
		}
		for (i = 0; i < params->nAttachedNCs; i++) {
			printf("Core %d is a %s-%s %s\n", i, params->attachedNC[i].ncAtomSymbol1, 
					params->attachedNC[i].ncAtomSymbol2, params->attachedNC[i].ncType);
			printf("The origin of this core is = %.2f % .2f % .2f\n", params->attachedNC[i].ncCenter.x,
						params->attachedNC[i].ncCenter.y, params->attachedNC[i].ncCenter.z);
			printf("The size of this core is = ");
			if (! strcmp(params->attachedNC[i].ncType, "QD")) {
				printf("D = %.2f A\n", params->attachedNC[i].ncSize.mag);
			}
			else if (! strcmp(params->attachedNC[i].ncType, "NR")) {
				printf("D = %.2f A ", params->attachedNC[i].ncSize.mag);
				printf("L = %.2f A\n", params->attachedNC[i].ncSize.z);
			}
			else if (! strcmp(params->attachedNC[i].ncType, "NPL")) {
				printf("L = %.2f A ", params->attachedNC[i].ncSize.x);
				printf("W = %.2f A ", params->attachedNC[i].ncSize.y);
				printf("H = %.2f A\n", params->attachedNC[i].ncSize.z);
			}
		}
		if (params->nStackingFaults) {
			printf("nStackingFaults = %d\n", params->nStackingFaults);
			for (j = 0; j < params->nStackingFaults; j++) {
				printf("stackFault[%d].type = %d\n", j, params->stackFault[j].type);
				printf("stackFault[%d].position = %d\n", j, params->stackFault[j].position);
			}
		}
	} 
	else {
		printf("confFile = %s\n", params->confFile);
		printf("confFileUnits = %s\n", params->confFileUnits);
		printf("buildNewNC = %d\n", params->buildNewNC);
		printf("ncType = %s\n", params->ncType);
		printf("nLayersToCut = %d\n", params->nLayersToCut);
		printf("remDanglingAtoms = %d\n", params->remDanglingAtoms);
	}
	printf("nMaxBonds = %d\n", params->nMaxBonds);
	printf("nMinBonds = %d\n", params->nMinBonds);
	printf("passivate = %d\n", params->passivate);
	writeSeparation(stdout);
	fflush(stdout);

	// Get the total number of atoms from first line of the configuration file
	if (! params->buildNewNC || ! params->buildNewBulk) {
		pf = fopen(params->confFile, "r");
		if (pf) {
			fscanf(pf, "%d", &params->nAtoms); 
			fclose(pf);
		}
		else {
			printf("PROGRAM EXITING: %s must be in current working directory\n", params->confFile);
			exit(EXIT_FAILURE);
		}
	}

	return;
}

/****************************************************************************/
// This function reads in the parameters of a nanocrystal structure

int fillNanocrystalStructure(nanostructure *nc, FILE *pf) {
	int i = 0;
	char field[1000], tmp[1000];

	while (fscanf(pf, "%s", field) != EOF) {
		if (! strcmp(field, "ncType")) fscanf(pf, "%s %s", tmp, &(nc->ncType));
		else if (! strcmp(field, "ncAtomSymbols")) fscanf(pf, "%s %s %s", tmp, &(nc->ncAtomSymbol1), &(nc->ncAtomSymbol2));
		else if (! strcmp(field, "ncCrystalStructure")) fscanf(pf, "%s %s", tmp, &(nc->ncCrystalStructure));
		else if (! strcmp(field, "ncCenter")) {
			fscanf(pf, "%s %lg %lg %lg", tmp, &(nc->ncCenter).x, &(nc->ncCenter).y, &(nc->ncCenter).z);
			nc->ncCenter.mag = retVectorMagnitude(nc->ncCenter);
		}
		else if (! strcmp(field, "ncSize")) {
			if (! strcmp(nc->ncType, "QD") || ! strcmp(nc->ncType, "attachedNCs")) { // qd diameter
				fscanf(pf, "%s %lg", tmp, &(nc->ncSize.mag)); 
			}
			else if (! strcmp(nc->ncType, "NR")) { // nr diameter and length
				fscanf(pf, "%s %lg %lg", tmp, &(nc->ncSize.mag), &(nc->ncSize.z)); 
			}
			else if (! strcmp(nc->ncType, "NPL")) { // npl length width and thickness
				fscanf(pf, "%s %lg %lg %lg", tmp, &(nc->ncSize.x), &(nc->ncSize.y), &(nc->ncSize.z));
			}
			else {
				printf("ncSize must be defined after ncType\n");
				printf("ncType = QD (NR, NPL or attachedNCs also allowed)\n");
				fflush(stdout);
				exit(EXIT_FAILURE);
			} 
		}
		else {
			printf("field = %s\n", field);
			printf("Invalid input field and/ or format following buildNewNC, nCores, nAttachedNCs\n");
			printf("Next lines must be (in any order, only one ncSize allowed):\n");
			printf("ncType = QD (NR, NPL or attachedNCs also allowed)\n");
			printf("ncAtomSymbols = Cd Se (Cd, Zn, Se, S, Te are allowed - Cd/Zn must come first)\n");
			printf("ncCrystalStructure = wurtzite (zincblende also allowed)\n");
			printf("ncCenter = xPos yPos zPos (in Angstroms)\n");
			printf("ncSize = qdDiameter (in Angstroms for ncType = QD or attachedNCs)\n");
			printf("ncSize = nrDiameter nrLength (in Angstroms for ncType = NR)\n");
			printf("ncSize = nplLength nplWidth nplThickness (in Angstroms for ncType = NPL)\n");
			fflush(stdout);
			exit(EXIT_FAILURE);
		}
		i++;
		if (i == 5) {
			return 0;
		}
	}

	return 1; // should not reach here
}


/****************************************************************************/

void readConf(atom *atoms, param *params) {
  FILE *pf;
  int i;

  params->nAtomTypes = 0; // initial value

  pf = fopen(params->confFile, "r"); 
  if (pf) {
    fscanf(pf, "%d", &i);        // skip the # of atoms (the first) line
    for (i = 0; i < params->nAtoms; i++) {
      fscanf(pf, "%s", atoms[i].symbol);
      if (isNewAtomType(atoms, i)) {
        atoms[i].type = params->nAtomTypes; // gives new integer for each new atom type
        params->nAtomTypes++;
      }
      fscanf(pf, "%lg %lg %lg", &atoms[i].pos.x, &atoms[i].pos.y, &atoms[i].pos.z);
      if (! strcmp(params->confFileUnits, "atomic")) atoms[i].pos = retScaledVector(atoms[i].pos, AUTOANG); // converted positions from a.u. to Ang
    }
    fclose(pf);
  }
  else {
    printf("PROGRAM EXITING: %s must be in current working directory\n", params->confFile);
    exit(EXIT_FAILURE);
  }

  return;
}

/****************************************************************************/
// Assigns all atoms in atoms a type and returns the number of unique atom types 

int assignAtomTypes(atom *atoms, int nAtoms) {
	int i, nAtomTypes = 0;

	for (i = 0; i < nAtoms; i++) atoms[i].type = 0;
    for (i = 0; i < nAtoms; i++) if (isNewAtomType(atoms, i)) {
        atoms[i].type = nAtomTypes; // gives new integer for each new atom type
        nAtomTypes++;
    }

	return nAtomTypes;
}

/****************************************************************************/
// returns 0/ false if atom symbol is already in the list or 1/true if 
// it is a new symbol in atoms[currIndex].symbol

int isNewAtomType(atom *atoms, int currIndex) {
  int i;

  for (i = 0; i < currIndex; i++) if (! strcmp(atoms[i].symbol, atoms[currIndex].symbol)) {
    atoms[currIndex].type = atoms[i].type; // make types equal since same symbols
    return 0; // not a new atom type
  }

  return 1; // new atom type
}

/****************************************************************************/
// Returns the atoms mass

double retAtomMass(char *atomSymbol) {

	if (! strcmp(atomSymbol, "Cd") || ! strcmp(atomSymbol, "Cdz")) return 112.411;
	else if (! strcmp(atomSymbol, "Zn") || ! strcmp(atomSymbol, "Znz")) return 65.380;
	else if (! strcmp(atomSymbol, "Se") || ! strcmp(atomSymbol, "Sez") || ! strcmp(atomSymbol, "Se1")) return 78.960;
	else if (! strcmp(atomSymbol, "S") || ! strcmp(atomSymbol, "Sz")) return 32.065;
	else if (! strcmp(atomSymbol, "Te")) return 127.60;
	else if (! strcmp(atomSymbol, "P1") || ! strcmp(atomSymbol, "P2")) return 1.0;
    else if (! strcmp(atomSymbol, "P3") || ! strcmp(atomSymbol, "P4")) return 1.0;
	else {
		writeSeparation(stdout);
		printf("WARNING: An atom symbol was found for which no atom mass is defined in retAtomMass\n");
		printf("The allowed atomic symbol are: Cd, Zn, Se, S, Te and P(1,2,3,4)\n");
		printf("The atomic symbol with no atomic mass is %s!\n", atomSymbol);
		writeSeparation(stdout);
		return 0.0;
	}

	return 0.0;
}

/****************************************************************************/
