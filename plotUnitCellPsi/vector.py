import numpy as np

def vecSum(r0, r1): 
	x_ans = r0[0] + r1[0]
	y_ans = r0[1] + r1[1]
	z_ans = r0[2] + r1[2]
	return np.array([x_ans, y_ans, z_ans])

def vecDiff(r0, r1): 
	x_ans = r0[0] - r1[0]
	y_ans = r0[1] - r1[1]
	z_ans = r0[2] - r1[2]
	return np.array([x_ans, y_ans, z_ans])

def vecDotProd(r0, r1): 
	x_ans = r0[0] - r1[0]
	y_ans = r0[1] - r1[1]
	z_ans = r0[2] - r1[2]
	return r0[0]*r1[0] + r0[1]*r1[1] + r0[2]*r1[2]

def vecLength(r):
	return np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)

def distance(r0, r1): 
	x_diff = r0[0] - r1[0]
	y_diff = r0[1] - r1[1]
	z_diff = r0[2] - r1[2]
	return np.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

def parallelWithTol(r0, r1, epsilon): 
	cos_angle = vecDotProd(r0, r1) / (vecLength(r0) * vecLength(r1))
	return (cos_angle<=1+epsilon) and (cos_angle>=1-epsilon)