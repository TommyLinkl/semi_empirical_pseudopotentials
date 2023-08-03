import os
import numpy as np
import random as rand

ANGTOAU = 1.88973
EPSILON = 0.000001

from vector import *

def read_xyz(filename):
	f = open(filename, "r")
	content = f.readlines()
	tot_num = int(content[0].strip())
	atomNames = np.array([content[line_num].strip().split()[0] for line_num in range(len(content)) if line_num>1])
	x_coord = np.array([float(content[line_num].strip().split()[1]) for line_num in range(len(content)) if line_num>1])
	y_coord = np.array([float(content[line_num].strip().split()[2]) for line_num in range(len(content)) if line_num>1])
	z_coord = np.array([float(content[line_num].strip().split()[3]) for line_num in range(len(content)) if line_num>1])
	coord = np.stack((x_coord, y_coord, z_coord), axis=-1)
	f.close()

	if (tot_num==len(atomNames)) and (tot_num==len(coord)):
		return (atomNames, coord)
	else:
		return ValueError("ERROR in reading .xyz file. Number of atoms not matching. ")

def read_cube_file(filename): 
	f = open(filename, "r")
	allLines = f.readlines()
	for i in range(len(allLines)): 
		allLines[i] = allLines[i].strip()

	# Header
	nAtoms = int(allLines[2].split()[0])
	origin = np.stack((float(allLines[2].split()[1]), float(allLines[2].split()[2]), float(allLines[2].split()[3])), axis=-1)
	nX = int(allLines[3].split()[0])
	vectorX = np.stack((float(allLines[3].split()[1]), float(allLines[3].split()[2]), float(allLines[3].split()[3])), axis=-1)
	nY = int(allLines[4].split()[0])
	vectorY = np.stack((float(allLines[4].split()[1]), float(allLines[4].split()[2]), float(allLines[4].split()[3])), axis=-1)
	nZ = int(allLines[5].split()[0])
	vectorZ = np.stack((float(allLines[5].split()[1]), float(allLines[5].split()[2]), float(allLines[5].split()[3])), axis=-1)

	# nAtoms, position of origin
	atom_type = np.array([int(allLines[line_num].split()[0]) for line_num in range(6, nAtoms+6)])
	atom_charge = np.array([float(allLines[line_num].split()[1]) for line_num in range(6, nAtoms+6)])
	x_coord = np.array([float(allLines[line_num].split()[2]) for line_num in range(6, nAtoms+6)])
	y_coord = np.array([float(allLines[line_num].split()[3]) for line_num in range(6, nAtoms+6)])
	z_coord = np.array([float(allLines[line_num].split()[4]) for line_num in range(6, nAtoms+6)])
	atom_coord = np.stack((x_coord, y_coord, z_coord), axis=-1)

	# volumetric_data
	volumetric_data = np.array(" ".join(allLines[nAtoms+6:]).strip().split()).astype(np.float)
	print(len(volumetric_data))
	print(int(allLines[3].split()[0]) * int(allLines[4].split()[0]) * int(allLines[5].split()[0]))

	volumetric_data = np.reshape(volumetric_data, (int(allLines[3].split()[0]), int(allLines[4].split()[0]), int(allLines[5].split()[0])))

	f.close()
	return (allLines[0:6], nAtoms, origin, nX, nY, nZ, vectorX, vectorY, vectorZ, atom_type, atom_charge, atom_coord, volumetric_data)

def write_UnitCell_cube_file(filename, firstSixLines, nAtoms, origin, nX, nY, nZ, vectorX, vectorY, vectorZ, atom_type, atom_charge, atom_coord, volumetric_data, xmin, ymin, zmin, xmax, ymax, zmax, buffer): 
	f = open(filename, "w")

	newNAtoms = 0
	for iAtom in range(nAtoms):
		if (atom_coord[iAtom][0]>=xmin-buffer) and (atom_coord[iAtom][0]<=xmax+buffer) and (atom_coord[iAtom][1]>=ymin-buffer) and (atom_coord[iAtom][1]<=ymax+buffer) and (atom_coord[iAtom][2]>=zmin-buffer) and (atom_coord[iAtom][2]<=zmax+buffer): 
			newNAtoms += 1

	if len(firstSixLines)!=6: 
		return ValueError("Header length is %d, incorrect. Should be 6." % len(firstSixLines))
	else: 
		f.write(firstSixLines[0])
		f.write("\n")
		f.write(firstSixLines[1])
		f.write("\n")
		f.write("%d   %.6f  %.6f  %.6f \n" % (newNAtoms, origin[0], origin[1], origin[2]))
		f.write(firstSixLines[3])
		f.write("\n")
		f.write(firstSixLines[4])
		f.write("\n")
		f.write(firstSixLines[5])
		f.write("\n")

	for iAtom in range(nAtoms):
		if (atom_coord[iAtom][0]>=xmin-buffer) and (atom_coord[iAtom][0]<=xmax+buffer) and (atom_coord[iAtom][1]>=ymin-buffer) and (atom_coord[iAtom][1]<=ymax+buffer) and (atom_coord[iAtom][2]>=zmin-buffer) and (atom_coord[iAtom][2]<=zmax+buffer): 
			f.write("%d  %.6f  %.6f  %.6f  %.6f \n" % (atom_type[iAtom], atom_charge[iAtom], atom_coord[iAtom][0], atom_coord[iAtom][1], atom_coord[iAtom][2]))

	for iX in range(nX):
		for iY in range(nY):
			for iZ in range(nZ): 
				vectorPoint = iX*vectorX+iY*vectorY+iZ*vectorZ + origin
				if (vectorPoint[0]>=xmin) and (vectorPoint[0]<=xmax) and (vectorPoint[1]>=ymin) and (vectorPoint[1]<=ymax) and (vectorPoint[2]>=zmin) and (vectorPoint[2]<=zmax):
					f.write("%.6f   " % volumetric_data[iX][iY][iZ])
				else: 
					f.write("0.000000   ")
				if (iZ%6==5):
					f.write("\n")
			f.write("\n")

	f.close()
	return 

def write_xyz(atomNames, coord, xyz_filename): 
	f = open(xyz_filename, "w")
	f.write("%d \n\n" % (len(atomNames)))
	for i in range(len(atomNames)):
		if (atomNames[i]=="P1") or (atomNames[i]=="P2") or (atomNames[i]=="P3"): 
			f.write("%s %f %f %f\n" % ("H", coord[i,0], coord[i,1], coord[i,2]))
		else:
			f.write("%s %f %f %f\n" % (atomNames[i], coord[i,0], coord[i,1], coord[i,2]))
	f.close()
	return

def write_unpassivated_xyz(atomNames, coord, xyz_filename): 
	f = open(xyz_filename, "w")
	f.write("%d \n # Comments \n" % (len(atomNames)-np.sum(atomNames=="P1")-np.sum(atomNames=="P2")-np.sum(atomNames=="P3")))
	for i in range(len(atomNames)):
		if (atomNames[i]!="P1") and (atomNames[i]!="P2") and (atomNames[i]!="P3"): 
			f.write("%s %f %f %f\n" % (atomNames[i], coord[i,0], coord[i,1], coord[i,2]))
	f.close()
	return

def write_dat(atomNames, coord, dat_filename): 
	f = open(dat_filename, "w")
	f.write("%d \n" % (len(atomNames)))
	for i in range(len(atomNames)):
		f.write("%s %f %f %f\n" % (atomNames[i], coord[i,0]*ANGTOAU, coord[i,1]*ANGTOAU, coord[i,2]*ANGTOAU))
	f.close()
	return

def retAtomMass(atomSymbol): 
	if atomSymbol=="In": 
		return 114.818
	elif atomSymbol=="Ga": 
		return 69.723
	elif atomSymbol=="As": 
		return 74.922
	elif atomSymbol=="P": 
		return 30.974
	else: 
		return 0.0

def write_lammpsconf(atomNames, coord, lammpsconf_filename): 
	f = open(lammpsconf_filename, "w")
	f.write("LAMMPS configuration data file\n\n")
	f.write("  %d atoms\n\n" % (len(atomNames)))
	f.write("  %d atom types\n\n" % (len(np.unique(atomNames))))
	f.write("  %.6f  %.6f  xlo xhi\n" % (coord.min(axis=0)[0]-20, coord.max(axis=0)[0]+20))
	f.write("  %.6f  %.6f  ylo yhi\n" % (coord.min(axis=0)[1]-20, coord.max(axis=0)[1]+20))
	f.write("  %.6f  %.6f  zlo zhi\n\n" % (coord.min(axis=0)[2]-20, coord.max(axis=0)[2]+20))
	f.write(" Masses\n\n")
	# print(np.unique(atomNames))
	for i in range(len(np.unique(atomNames))): 
		for j in range(len(atomNames)): 
			if atomNames[j]==np.unique(atomNames)[i]: 
				f.write("  %d   %f\n" % (i+1, retAtomMass(atomNames[j])))
				break
	f.write("\n Atoms\n\n")
	for i in range(len(atomNames)): 
		for j in range(len(np.unique(atomNames))): 
			if atomNames[i]==np.unique(atomNames)[j]: 
				f.write("        %d    %d     %.8f  %.8f  %.8f\n" % (i+1, j+1, coord[i,0], coord[i,1], coord[i,2]))
				# print(atomNames[i], j+1)
	f.close()
	return
