# # # # Dipti Jasrasaria
# January 2022
# calculates reference tetrahedron volumes

import numpy as np

# tetrahedron with bond length of 1, side length of sqrt(8/3), center at origin
v1 = np.array([np.sqrt(8./9.), 0., -1./3.])
v2 = np.array([-np.sqrt(2./9.), np.sqrt(2./3.), -1./3.])
v3 = np.array([-np.sqrt(2./9.), -np.sqrt(2./3.), -1./3.])
v4 = np.array([0., 0., 1.])

def vol(a, b, c, d):
    return np.abs(np.dot(a-d, np.cross(b-d, c-d))/6.)

# # # # wurtzite # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cd_se = 2.6326 # Cd-Se bond length (Ang)
cd_s = 2.5292 # Cd-S bond length (Ang)

# Se center, 4 Cd vertices
vol_se4cd = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_se)

# S center, 4 Cd vertices
vol_s4cd = vol(v1*cd_s, v2*cd_s, v3*cd_s, v4*cd_s)

# Cd center, 4 Se vertices
vol_cd4se = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_se)

# Cd center, 3 Se vertices, 1 S vertex
vol_cd3se1s = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_s)

# Cd center, 2 Se vertices, 2 S vertices
vol_cd2se2s = vol(v1*cd_se, v2*cd_se, v3*cd_s, v4*cd_s)

# Cd center, 1 Se vertex, 3 S vertices
vol_cd1se3s = vol(v1*cd_se, v2*cd_s, v3*cd_s, v4*cd_s)

# Cd center, 4 S vertices
vol_cd4s = vol(v1*cd_s, v2*cd_s, v3*cd_s, v4*cd_s)

print()
# # # # zincblende # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cd_se = 2.6233 # Cd-Se bond length (Ang)
cd_s = 2.5193 # Cd-S bond length (Ang)
in_p = 2.5228

# Se center, 4 Cd vertices
vol_se4cd = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_se)

# S center, 4 Cd vertices
vol_s4cd = vol(v1*cd_s, v2*cd_s, v3*cd_s, v4*cd_s)

# Cd center, 4 Se vertices
vol_cd4se = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_se)

# Cd center, 3 Se vertices, 1 S vertex
vol_cd3se1s = vol(v1*cd_se, v2*cd_se, v3*cd_se, v4*cd_s)

# Cd center, 2 Se vertices, 2 S vertices
vol_cd2se2s = vol(v1*cd_se, v2*cd_se, v3*cd_s, v4*cd_s)

# Cd center, 1 Se vertex, 3 S vertices
vol_cd1se3s = vol(v1*cd_se, v2*cd_s, v3*cd_s, v4*cd_s)

# Cd center, 4 S vertices
vol_cd4s = vol(v1*cd_s, v2*cd_s, v3*cd_s, v4*cd_s)

vol_in4p = vol(v1*in_p, v2*in_p, v3*in_p, v4*in_p)
print(vol_in4p)

'''
# # # # verify derivatives # # # # # # # # # # # # # # # # # # # # # # # # # # #
dx = 0.013299653

dxx = np.array([dx, 0., 0.])
dxy = np.array([0., dx, 0.])
dxz = np.array([0., 0., dx])

if np.dot(v1-v4, np.cross(v2-v4, v3-v4)) < 0:
    n = -1.
else:
    n = 1.

# # # # dv / dv1
dvdv1x = (vol(v1+dxx,v2,v3,v4) - vol(v1-dxx,v2,v3,v4))/(2.*dx)
dvdv1y = (vol(v1+dxy,v2,v3,v4) - vol(v1-dxy,v2,v3,v4))/(2.*dx)
dvdv1z = (vol(v1+dxz,v2,v3,v4) - vol(v1-dxz,v2,v3,v4))/(2.*dx)
tmp = (n/6.)*np.cross(v2-v4, v3-v4)

print('dx: ')
print(dx)
print('dv / dv1x: ')
print(dvdv1x)
print(tmp[0])
print('dv / dv1y: ')
print(dvdv1y)
print(tmp[1])
print('dv / dv1z: ')
print(dvdv1z)
print(tmp[2])

# # # # dv / dv2
dvdv2x = (vol(v1,v2+dxx,v3,v4) - vol(v1,v2-dxx,v3,v4))/(2.*dx)
dvdv2y = (vol(v1,v2+dxy,v3,v4) - vol(v1,v2-dxy,v3,v4))/(2.*dx)
dvdv2z = (vol(v1,v2+dxz,v3,v4) - vol(v1,v2-dxz,v3,v4))/(2.*dx)

print('dv / dv2x: ')
print(dvdv2x)
print((n/6.)*np.dot(v1-v4, np.cross(np.array([1., 0., 0.]), v3-v4)))
print('dv / dv2y: ')
print(dvdv2y)
print((n/6.)*np.dot(v1-v4, np.cross(np.array([0., 1., 0.]), v3-v4)))
print('dv / dv2z: ')
print(dvdv2z)
print((n/6.)*np.dot(v1-v4, np.cross(np.array([0., 0., 1.]), v3-v4)))

# # # # dv / dv3
dvdv3x = (vol(v1,v2,v3+dxx,v4) - vol(v1,v2,v3-dxx,v4))/(2.*dx)
dvdv3y = (vol(v1,v2,v3+dxy,v4) - vol(v1,v2,v3-dxy,v4))/(2.*dx)
dvdv3z = (vol(v1,v2,v3+dxz,v4) - vol(v1,v2,v3-dxz,v4))/(2.*dx)

print('dv / dv3x: ')
print(dvdv3x)
print((n/6.)*np.dot(v1-v4, np.cross(v2-v4, np.array([1., 0., 0.]))))
print('dv / dv3y: ')
print(dvdv3y)
print((n/6.)*np.dot(v1-v4, np.cross(v2-v4, np.array([0., 1., 0.]))))
print('dv / dv3z: ')
print(dvdv3z)
print((n/6.)*np.dot(v1-v4, np.cross(v2-v4, np.array([0., 0., 1.]))))

# # # # dv / dv4
dvdv4x = (vol(v1,v2,v3,v4+dxx) - vol(v1,v2,v3,v4-dxx))/(2.*dx)
dvdv4y = (vol(v1,v2,v3,v4+dxy) - vol(v1,v2,v3,v4-dxy))/(2.*dx)
dvdv4z = (vol(v1,v2,v3,v4+dxz) - vol(v1,v2,v3,v4-dxz))/(2.*dx)

print('dv / dv4x: ')
print(dvdv4x)
print(-(dvdv1x+dvdv2x+dvdv3x))
print('dv / dv4y: ')
print(dvdv4y)
print(-(dvdv1y+dvdv2y+dvdv3y))
print('dv / dv4z: ')
print(dvdv4z)
print(-(dvdv1z+dvdv2z+dvdv3z))
'''