from vector import *
from truncate_unitCell import * 

(headers, nAtoms, origin, nX, nY, nZ, vectorX, vectorY, vectorZ, atom_type, atom_charge, atom_coord, volumetric_data) = read_cube_file("psiElec-ni-0.cube")

write_UnitCell_cube_file("UnitCell_psiElec-ni-0.cube", headers, nAtoms, origin, nX, nY, nZ, vectorX, vectorY, vectorZ, atom_type, atom_charge, atom_coord, volumetric_data, -5, -5, -5, 5, 5, 5, 0.7)