#######################################################################################################
# This python script creates plots from mco files.
#######################################################################################################

import sys
import os
from enum import Enum
import numpy as NUM
import matplotlib.pyplot as PLT
import argparse


########################################################################################

parser = argparse.ArgumentParser(description='Plot mco files')

parser.add_argument('files', type=str, nargs='+',
	help='list of mco files')

parser.add_argument('--compare', dest='COMPARE_MODE_ENABLED', action='store_const', const=True, default=False,
	help='Compare given mco files by using same plot for corresponding data')

parser.add_argument('--boundaries', dest='PLOT_BOUNDARIES', action='store_const', const=True, default=False,
	help='Plot boundaries on second y axis')

parser.add_argument('-o', dest='outFolder', action='store', default='comparison',
	help='Output folder')

args = parser.parse_args()


########################################################################################

# Helper type for parsing
# We can parse and plot the following quantities:
# - A_l: Absorption as function of layer [-] 
# - A_z: Absorption as function of depth [1/cm] 
# - Rd_r: Diffuse reflectance as function of radius [1/cm2]
# - Rd_a: Diffuse reflectance as function of exit angle [1/sr] 
# - Tt_r: Total transmittance (at last layer) as function of radius [1/cm2] 
# - Tt_a: Total transmittance as function of exit angle [1/sr] 
# - A_rz: 2D probability density in turbid media over r & z. [1/cm3]
# - Rd_ra: 2D distribution of diffuse reflectance. [1/(cm2 sr)]
# - Tt_ra: 2D distribution of total transmittance. [1/(cm2 sr)]
Plottable = Enum('Plottable', 'NONE A_l A_z Rd_r Rd_a Tt_r Tt_a A_rz Rd_ra Tt_ra')

########################################################################################


# Plot and label 2- or 3-dimensional data array
# Return PDF file path to be saved
def plotToPDF_internal(folder, components, data, resolutions, units, boundaries_x, boundaries_y):

	path = ''

	if len(components) == 2:
		if args.PLOT_BOUNDARIES:

			# Make 2 y-axes
			fig, ax1 = PLT.subplots()
			ax1.plot(data)
			ax1.set_xlabel(components[1]+' ['+str(resolutions[0])+units[0]+']')
			ax1.set_ylabel(components[0]+' ['+str(resolutions[1])+units[1]+']')
			ax1.set_yscale('log')
			ax2 = ax1.twinx()
			ax2.set_ylabel('height [cm]')
			for i in range(0, len(boundaries_x)):

				# If boundary doesnt fill the domain, continue it as straight line
				if boundaries_x[i][len(boundaries_x[i])-1]/resolutions[0] < len(data):
					boundaries_x[i].append(len(data)*resolutions[0])
					boundaries_y[i].append(boundaries_y[i][len(boundaries_y[i])-1])

				ax2.plot(list(map(lambda x: x/resolutions[0], boundaries_x[i])), boundaries_y[i], linestyle='dashed', linewidth=1)
		else:
			PLT.plot(data)
			PLT.ylabel(components[0]+' ['+str(resolutions[1])+units[1]+']')
			PLT.yscale('log')
			PLT.xlabel(components[1]+' ['+str(resolutions[0])+units[0]+']')
		PLT.title(components[0]+'_'+components[1])
		path = folder + '/' + components[0]+'_'+components[1] + '.pdf'
	
	elif len(components) == 3:
		PLT.imshow(NUM.matrix(data).transpose())
		PLT.title(components[0]+'_'+components[1]+'_'+components[2])
		PLT.xlabel(components[1]+' ['+str(resolutions[0])+units[0]+']')
		PLT.ylabel(components[2]+' ['+str(resolutions[1])+units[1]+']')
		#TODO label for 3rd dimension (unit of radiance)
		PLT.colorbar()
		path = folder + '/' + components[0]+'_'+components[1]+'_'+components[2] + '.pdf'

		#FIXME For now, cut away low pixels to see the few actual data points, 
		# but this is not good anymore when grid size changes.
		# Please fix this when actually using the 3D data.
		PLT.ylim(len(data[0])*0.1, 0)
		PLT.xlim(0, len(data)*0.01)

	else:
		print('Cannot plot ' + str(len(components)) + '-dimensional data')

	return path

# Plot one 2- or 3-dimensional data array
def plotToPDF(folder, components, data, resolutions, units, boundaries_x, boundaries_y):

	PLT.figure()

	path = plotToPDF_internal(folder, components, data, resolutions, units, boundaries_x, boundaries_y)

	if args.PLOT_BOUNDARIES:
		boundary_names = []
		for i in range(0, len(boundaries_x)):
			boundary_names.append('Boundary ' + str(i))
		PLT.legend(boundary_names)

	if not os.path.exists(folder): os.makedirs(folder)
	PLT.savefig(path)
	PLT.close()
	
	print('Written ' + path)

# Plot a batch of data arrays
def plotAllToPDF(folder, components, data_arrays, resolutions, units, data_names, boundaries_x, boundaries_y):

	PLT.figure()

	path = ''

	for data in data_arrays:
		path = plotToPDF_internal(folder, components, data, resolutions, units, boundaries_x, boundaries_y)

	boundary_names = []
	if args.PLOT_BOUNDARIES:
		for i in range(0, len(boundaries_x)):
			boundary_names.append('Boundary ' + str(i))
	PLT.legend(data_names + boundary_names)

	if not os.path.exists(folder): os.makedirs(folder)
	PLT.savefig(path)
	PLT.close()
	
	print('Written ' + path)


########################################################################################


#
def parseInParams(stripped_lines):

	n_layers = 0

	dz = 0
	dr = 0

	nz = 0
	nr = 0
	na = 0

	for i in range(0, len(stripped_lines)):
		l = stripped_lines[i]
		if l[0:3]=='InP':
			deltas = stripped_lines[i+3].split()
			dz = float(deltas[0])
			dr = float(deltas[1])
			dims = stripped_lines[i+4].split()
			nz = int(dims[0])
			nr = int(dims[1])
			na = int(dims[2])
			n_layers = int(stripped_lines[i+5].split()[0])
			break

	# Check
	if n_layers==0 or dz==0 or dr==0 or nz==0 or nr==0 or na==0:
		print('layer and detection arrays must have non-zero size ('+mco+')')
		exit()

	return n_layers, dz, dr, nz, nr, na


#
def parseOutData(stripped_lines):

	Rd_r = []
	Rd_a = []
	Rd_ra = [[]]

	A_l = []
	A_z = []
	A_rz = [[]]

	Tt_r = []
	Tt_a = []
	Tt_ra = [[]]

	current_plottable = Plottable.NONE

	for l in stripped_lines:

		if l[0:3]=='A_l':
			current_plottable = Plottable.A_l
			continue
		elif l[0:3]=='A_z':
			current_plottable = Plottable.A_z
			continue
		elif l[0:4]=='A_rz':
			current_plottable = Plottable.A_rz
			continue
		elif l[0:5]=='Rd_ra':
			current_plottable = Plottable.Rd_ra
			continue
		elif l[0:5]=='Tt_ra':
			current_plottable = Plottable.Tt_ra
			continue
		elif l[0:4]=='Rd_r':
			current_plottable = Plottable.Rd_r
			continue
		elif l[0:4]=='Rd_a':
			current_plottable = Plottable.Rd_a
			continue
		elif l[0:4]=='Tt_r':
			current_plottable = Plottable.Tt_r
			continue
		elif l[0:4]=='Tt_a':
			current_plottable = Plottable.Tt_a
			continue

		if current_plottable==Plottable.A_l:
			A_l.append(float(l))
		elif current_plottable==Plottable.A_z:
			A_z.append(float(l))
		elif current_plottable==Plottable.Rd_r:
			Rd_r.append(float(l))
		elif current_plottable==Plottable.Rd_a:
			Rd_a.append(float(l))
		elif current_plottable==Plottable.Tt_r:
			Tt_r.append(float(l))
		elif current_plottable==Plottable.Tt_a:
			Tt_a.append(float(l))
		elif current_plottable==Plottable.A_rz:
			numbers = l.split()
			for i in range(0, len(numbers)):
				if len(A_rz[len(A_rz)-1]) < nz:
					A_rz[len(A_rz)-1].append( float(numbers[i]) )
				else:
					A_rz.append([ float(numbers[i]) ])
		elif current_plottable==Plottable.Rd_ra:
			numbers = l.split()
			for i in range(0, len(numbers)):
				if len(Rd_ra[len(Rd_ra)-1]) < na:
					Rd_ra[len(Rd_ra)-1].append( float(numbers[i]) )
				else:
					Rd_ra.append([ float(numbers[i]) ])
		elif current_plottable==Plottable.Tt_ra:
			numbers = l.split()
			for i in range(0, len(numbers)):
				if len(Tt_ra[len(Tt_ra)-1]) < na:
					Tt_ra[len(Tt_ra)-1].append( float(numbers[i]) )
				else:
					Tt_ra.append([ float(numbers[i]) ])

	return Rd_r, Rd_a, Rd_ra, A_l, A_z, A_rz, Tt_r, Tt_a, Tt_ra

#
def checkOutData(Rd_r, Rd_a, Rd_ra, A_l, A_z, A_rz, Tt_r, Tt_a, Tt_ra, nr, nz, na):

	if len(Rd_r) != nr:
		print('Incorrect data dimensions in Rd_r')
		exit()
	if len(Rd_ra) != nr:
		print('Incorrect data dimensions in Rd_ra')
		exit()
	for l in Rd_ra:
		if len(l) != na:
			print('Incorrect data dimensions in Rd_ra')
			exit()

	if len(A_rz) != nr:
		print('Incorrect data dimensions in A_rz')
		exit()
	for l in A_rz:
		if len(l) != nz:
			print('Incorrect data dimensions in A_rz')
			exit()

	if len(Tt_ra) != nr:
		print('Incorrect data dimensions in Tt_ra')
		exit()
	for l in Tt_ra:
		if len(l) != na:
			print('Incorrect data dimensions in Tt_ra')
			exit()


def parseBoundaries(stripped_lines):

	all_x = []
	all_y = []

	for l in stripped_lines:
		if l[0] == 'b':
			b = l.split()
			b.pop(0)
			n = b.pop(0)

			x = []
			y = []

			current_x = 0

			for i in range(0, len(b)):
				if i%2==0:
					y.append(float(b[i]))
				else:
					current_x += float(b[i])
					x.append(current_x)

			all_x.append(x)
			all_y.append(y)

	return all_x, all_y

########################################################################################

# Multiplies distance in x with values in array
def multiplyXOntoY(array, step):
	result = []
	d = 0
	for v in array:
		result.append(v * d)
		d += step
	return result

########################################################################################

mcos = []

Rd_r_arrays = []
Rd_a_arrays = []
Tt_r_arrays = []
Tt_a_arrays = []
last_na = 0
last_nr = 0
last_dr = 0

last_boundaries_x = []
last_boundaries_y = []
all_boundaries_x = {}
all_boundaries_y = {}

# Extract the data from each mco file
for i in range(0, len(args.files)):

	mco = args.files[i]

	print('Parsing ' + mco + '...')

	with open(mco) as f:

		lines = f.read().splitlines()
		stripped_lines = [l.strip() for l in lines if l.strip()!='' and l[0]!='#']


		#TODO handle IGNORE_A case


		n_layers, dz, dr, nz, nr, na = parseInParams(stripped_lines)

		Rd_r, Rd_a, Rd_ra, A_l, A_z, A_rz, Tt_r, Tt_a, Tt_ra = parseOutData(stripped_lines)

		checkOutData(Rd_r, Rd_a, Rd_ra, A_l, A_z, A_rz, Tt_r, Tt_a, Tt_ra, nr, nz, na)

		boundaries_x, boundaries_y = parseBoundaries(stripped_lines)
		all_boundaries_x[mco] = boundaries_x
		all_boundaries_y[mco] = boundaries_y

		if args.COMPARE_MODE_ENABLED:

			mcos.append(mco)

			# For now, only interested in 2-dimensional R and T data

			Rd_r_arrays.append(Rd_r)
			Rd_a_arrays.append(Rd_a)
			Tt_r_arrays.append(Tt_r)
			Tt_a_arrays.append(Tt_a)

			if i > 0 and (na != last_na or nr != last_nr or dr != last_dr):
				print('Cannot compare data with different scales')
				exit()

			last_na = na
			last_nr = nr
			last_dr = dr

			if args.PLOT_BOUNDARIES:
				if i > 0:
					if len(last_boundaries_x) != len(boundaries_x) or len(last_boundaries_y) != len(boundaries_y):
						print('Cannot plot different number of boundaries in compare mode')
						exit()
					for j in range(0, len(boundaries_x)):
						if len(last_boundaries_x[j]) != len(boundaries_x[j]):
							print('Cannot plot boundaries with different number of samples in compare mode')
							exit()
						for k in range(0, len(boundaries_x[j])):
							if last_boundaries_x[j][k] != boundaries_x[j][k]:
								print('Cannot plot different boundaries in compare mode')
								exit()
					for j in range(0, len(boundaries_y)):
						if len(last_boundaries_x[j]) != len(boundaries_y[j]):
							print('Cannot plot boundaries with different number of samples in compare mode')
							exit()
						for k in range(0, len(boundaries_y[j])):
							if last_boundaries_x[j][k] != boundaries_y[j][k]:
								print('Cannot plot different boundaries in compare mode')
								exit()
				last_boundaries_x = boundaries_x
				last_boundaries_y = boundaries_y

		else:

			folder = mco[0:-4]
			plotToPDF(folder, ['Rd_times_r', 'r'], multiplyXOntoY(Rd_r, dr), [dr, 1],          ['cm', 'cm-2 * cm'], # read: radius in cm and reflectance in J per square cm multiplied with r
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Rd', 'r'],         Rd_r,                     [dr, 1],          ['cm', 'cm-2'], # read: radius in cm and reflectance in J per square cm multiplied with r
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Rd', 'a'],         Rd_a,                     [na/90.0, 1],     ['deg', 'sr-1'], # reflectance per solid angle
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Rd', 'r', 'a'],    Rd_ra,                    [dr, na/90.0, 1], ['cm', 'deg', 'cm-2 * sr-1'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['A', 'l'],          A_l,                      [1, 1],           ['layer', '-'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['A', 'z'],          A_z,                      [dz, 1],          ['cm', 'cm-1'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['A', 'r', 'z'],     A_rz,                     [dr, dz, 1],      ['cm', 'cm', 'cm-3'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Tt', 'r', 'a'],    Tt_ra,                    [dr, na/90.0, 1], ['cm', 'deg', 'cm-2 * sr-1'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Tt_times_r', 'r'], multiplyXOntoY(Tt_r, dr), [dr, 1],          ['cm', 'cm-2 * cm'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Tt', 'r'],         Tt_r,                     [dr, 1],          ['cm', 'cm-2'],
				boundaries_x, boundaries_y)
			plotToPDF(folder, ['Tt', 'a'],         Tt_a,                     [na/90.0, 1],     ['deg', 'sr-1'],
				boundaries_x, boundaries_y)


if args.COMPARE_MODE_ENABLED:

	mcos = list(map(lambda mco: os.path.split(mco)[1], mcos))

	plotAllToPDF(args.outFolder, ['Rd_times_r', 'r'], map(lambda array: multiplyXOntoY(array, last_dr), Rd_r_arrays), [last_dr, 1],      ['cm', 'cm-2 * cm'], mcos, boundaries_x, boundaries_y)
	plotAllToPDF(args.outFolder, ['Rd', 'r'],         Rd_r_arrays,                                                    [last_dr, 1],      ['cm', 'cm-2'], mcos, boundaries_x, boundaries_y)
	plotAllToPDF(args.outFolder, ['Rd', 'a'],         Rd_a_arrays,                                                    [last_na/90.0, 1], ['deg', 'sr-1'], mcos, boundaries_x, boundaries_y)
	plotAllToPDF(args.outFolder, ['Tt_times_r', 'r'], map(lambda array: multiplyXOntoY(array, last_dr), Tt_r_arrays), [last_dr, 1],      ['cm', 'cm-2 * cm'], mcos, boundaries_x, boundaries_y)
	plotAllToPDF(args.outFolder, ['Tt', 'r'],         Tt_r_arrays,                                                    [last_dr, 1],      ['cm', 'cm-2'], mcos, boundaries_x, boundaries_y)
	plotAllToPDF(args.outFolder, ['Tt', 'a'],         Tt_a_arrays,                                                    [last_na/90.0, 1], ['deg', 'sr-1'], mcos, boundaries_x, boundaries_y)