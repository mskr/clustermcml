#######################################################################################################
# This python script generates 1D representations of radial heightfields, to use as MCML boundaries.
#######################################################################################################

import numpy as NUM
import matplotlib.pyplot as PLT
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


# Uncomment to test generators visually
# PLT.figure(1)
# PLT.plot( generateCurve(-5, 50, 0.01)[0::2] ) # plot only every 2nd element, i.e. height values
# PLT.show()