#######################################################################################################
# This python script generates 1D representations of radial heightfields, to use as MCML boundaries.
#######################################################################################################

import numpy as NUM
import matplotlib.pyplot as PLT
import math


def centerAroundMean(heightfield):
	nsamples = int(len(heightfield)/2)
	h_sum = 0
	for i in range(0, nsamples):
		h_sum = h_sum + heightfield[i*2]
	h_mean = h_sum / nsamples
	for i in range(0, nsamples):
		heightfield[i*2] = heightfield[i*2] - h_mean
	return heightfield


def generateZeros(nsamples, spacing):
	heightfield = []
	for i in range(0, nsamples):
		heightfield.append(0)
		heightfield.append(spacing)
	return heightfield


# Generate curve as small cut of the upper part of a circle,
# roughly looking like this:
#   __
#     ''''----____
#
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
# PLT.plot( generateCurve(-2, 50, 0.01)[0::2] ) # plot only every 2nd element, i.e. height values
# PLT.show()