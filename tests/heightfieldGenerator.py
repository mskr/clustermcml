#######################################################################################################
# This python script generates 1D representations of radial heightfields, to use as MCML boundaries.
#######################################################################################################

import numpy as NUM
import matplotlib.pyplot as PLT
from mpl_toolkits.mplot3d import axes3d
import math


# Generate flat heightfield of given size with given uniform spacing
def generateZeros(nsamples, spacing):
	heightfield = []
	for i in range(0, nsamples):
		heightfield.append(0)
		heightfield.append(spacing)
	return heightfield

# Center heightfield to be "comparable" with flat surface
def centerAroundMean(heightfield):
	nsamples = int(len(heightfield)/2)

	# Get mean height, with respect to area.
	# Note that outer values contribute to a larger part of the radial surface.
	r = 0
	numer = 0
	denom = 0
	for i in range(0, nsamples):
		r += heightfield[i*2+1]
		numer += heightfield[i*2] * 2*math.pi*r
		denom += 2*math.pi*r

	h_mean = numer / denom

	for i in range(0, nsamples):
		heightfield[i*2] -= h_mean

	return heightfield


# Generate curve as small cut of the upper part of a circle,
# roughly looking like this:
#   __
#     ''''----____
#
# Spacing is uniform.
def generateCurve(maxheight, nsamples, spacing):
	a_max = math.pi/4 # maximum angle (end of curve)

	y_base = math.cos(a_max)

	norm = 1 / (1 - y_base)

	heightfield = []

	a_step = a_max / nsamples

	for i in range(1, nsamples+1):
		y = math.cos(i * a_step) - y_base
		y = y * norm * maxheight
		heightfield.append(y)
		heightfield.append(spacing)

	return centerAroundMean(heightfield)


# Generate uniform noise with uniform spacing
def generateNoise(maxheight, nsamples, spacing):
	heights = NUM.random.sample(nsamples).tolist()
	heights = list(map(lambda h : h * maxheight, heights))
	spacings = [spacing] * nsamples
	heightfield = []
	for i in range(0, nsamples):
		heightfield.append(heights[i])
		heightfield.append(spacings[i])
	return centerAroundMean(heightfield)


###########################################################################################################


def lerp(v0, v1, a):
	return v0 + (v1 - v0) * a

def getHeightAt(radial_heightfield, pos):
	heights = radial_heightfield[0::2]
	spacings = radial_heightfield[1::2]
	n = len(heights)
	i = 0
	r = 0.0
	r_last = 0.0
	while r < pos:
		r_last = r;
		r += spacings[i if i < n else n-1]
		i+=1
		if i >= n: return heights[n-1] # pos outside sampled range, return last sample
	# Now pos is between r_last and r
	a = 0.0 if r==r_last else (pos - r_last) / (r - r_last)
	h0 = heights[0 if i==0 else i-1 if i<=n else n-1] # read sample, clamping to first resp. last sample
	h1 = heights[i if i<n else n-1] # read sample, clamping to last sample
	return lerp(h0, h1, a)

def getZ(radial_heightfield, X, Y):
	Z = []
	for i in range(0, len(X)):
		Z.append([])
		for j in range(0, len(X[i])):
			Z[i].append(getHeightAt(radial_heightfield, NUM.linalg.norm([X[i][j], Y[i][j]])))
	return NUM.asarray(Z)

def get3DMesh(radial_heightfield):
	heights = radial_heightfield[0::2]
	spacings = radial_heightfield[1::2]
	radii = []
	r = 0
	for s in spacings:
		radii.append(r)
		r += s

	x = y = NUM.arange(-r-.2, r+.2, min(spacings))
	X, Y = NUM.meshgrid(x, y)

	return X, Y, getZ(radial_heightfield, X, Y)


def get2DPoints(radial_heightfield):
	y = radial_heightfield[0::2]
	x = []
	r = 0
	for s in radial_heightfield[1::2]:
		x.append(r)
		r += s
	return x, y


###########################################################################################################


def plot2D():
	PLT.figure(1)
	PLT.plot( generateCurve(-5, 8, 0.1)[0::2] ) # plot only every 2nd element, i.e. height values
	PLT.show()


def plot3D():

	hfield = generateCurve(1, 10, 0.1)

	X, Y, Z = get3DMesh(hfield)

	fig = PLT.figure(figsize=(8,4))
	ax = fig.add_subplot(122, projection='3d')

	ax.plot_wireframe(X, Y, Z)

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	x, y = get2DPoints(hfield)
	ax2 = fig.add_subplot(121, xticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
	ax2.scatter(x, y)
	ax2.set_xlabel('r')
	ax2.set_ylabel('Z')

	PLT.show()#savefig('heightfield.pdf')

# plot3D() # uncomment to look at heightfield plots


###########################################################################################################


# curve = generateCurve(-5, 8, 0.1)
# heights = curve[0::2]
# spacings = curve[1::2]

# radii = []
# r = 0
# for s in spacings:
# 	radii.append(r)
# 	r += s

# X = []
# Y = []
# Z = []
# for i in reversed(range(0,len(radii))):
# 	r = radii[i]
# 	last_r = radii[i-1] if i >= 0 else 0
# 	h = heights[i]
# 	last_h = heights[i-1] if i >= 0 else heights[i]
# 	a = 0
# 	while a < 2*math.pi:

# 		next_a = a + 2*math.pi/16

# 		X.append([r*math.sin(a), r*math.sin(next_a)])
# 		Y.append([r*math.cos(a), r*math.cos(next_a)])
# 		Z.append([h,h])

# 		X.append([r*math.sin(next_a), last_r*math.sin(next_a)])
# 		if next_a < 2*math.pi:
# 			X.append([last_r*math.sin(next_a), r*math.sin(next_a)])

# 		Y.append([r*math.cos(next_a), last_r*math.cos(next_a)])
# 		if next_a < 2*math.pi:
# 			Y.append([last_r*math.cos(next_a), r*math.cos(next_a)])
		
# 		Z.append([h,last_h])
# 		if next_a < 2*math.pi:
# 			Z.append([last_h,h])

# 		a = next_a

# X = NUM.asarray(X)
# Y = NUM.asarray(Y)
# Z = NUM.asarray(Z)
