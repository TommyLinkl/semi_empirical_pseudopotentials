import numpy as np
from numpy.linalg import svd

'''
a = np.random.randn(9, 6)
print(a)

u, s, vh = np.linalg.svd(a, full_matrices=True)
print(u.shape, s.shape, vh.shape)
print(u)
print(s)
print(vh)

u, s, vh = np.linalg.svd(a, full_matrices=False)
print(u.shape, s.shape, vh.shape)
print(u)
print(s)
print(vh)
'''

# read

# tmp = np.array([float(line.strip().split()[0]) for line in open("DEBUG_psitot.dat", 'r')])

f = open("DEBUG_psitot.dat", "r")
nline = 0
while True:
	line = f.readline()
	if not line:
		break
	else: 
		# print(nline)
		tmp = np.array([float(x) for x in line.split()])
		# print(tmp.shape)
		if (nline==0): 
			matrix = np.copy(tmp)
		else: 
			matrix = np.vstack((matrix, tmp))
		nline += 1
		# print(matrix.shape)
f.close()

matrix = np.transpose(matrix)

# svd
try: 
	u, s, vh = np.linalg.svd(matrix, full_matrices=False)
except LinAlgError:
    print("LinAlgError")

print(u.shape, s.shape, vh.shape)
# print(u)
print(s)
# print(vh)


test_matrix = np.array([[1, 0, 4], [2, 1, 0], [3,1,5]])
print(test_matrix)
# svd
try: 
	u, s, vh = np.linalg.svd(test_matrix, full_matrices=False)
except LinAlgError:
    print("LinAlgError")

print(u.shape, s.shape, vh.shape)
# print(u)
print(s)
# print(vh)
