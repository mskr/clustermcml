#######################################################################################################
# This python script creates plots from mco files.
#######################################################################################################

import sys
import os
from enum import Enum
import numpy as NUM
import matplotlib.pyplot as PLT

if len(sys.argv) < 2:
	print('Usage: python plotter.py <list of mco files>')
	exit()

Plottable = Enum('Plottable', 'NONE A_l A_z Rd_r Rd_a Tt_r Tt_a A_rz Rd_ra Tt_ra')

# This function can plot and save data to a file...
def plotToPDF(folder, components, data, resolutions, units):

	PLT.figure()

	if len(components) == 2:
		PLT.plot(data)
		PLT.title(components[0]+'_'+components[1])
		PLT.xlabel(components[1]+' ['+str(resolutions[0])+units[0]+']')
		PLT.ylabel(components[0])
		PLT.yscale('log')
		path = folder + '/' + components[0]+'_'+components[1] + '.pdf'
	
	elif len(components) == 3:
		PLT.imshow(NUM.matrix(data).transpose())
		PLT.title(components[0]+'_'+components[1]+'_'+components[2])
		PLT.xlabel(components[1]+' ['+str(resolutions[0])+units[0]+']')
		PLT.ylabel(components[2]+' ['+str(resolutions[1])+units[1]+']')
		PLT.colorbar()
		path = folder + '/' + components[0]+'_'+components[1]+'_'+components[2] + '.pdf'

		# cut away low pixels to see the few actual data points
		PLT.ylim(len(data[0])*0.1, 0)
		PLT.xlim(0, len(data)*0.01)

	else:
		print('Cannot plot ' + str(len(components)) + '-dimensional data')

	if not os.path.exists(folder): os.makedirs(folder)
	PLT.savefig(path)
	PLT.close()
	
	print('Written ' + path)


# ... but first we need to extract the data from each mco file
for i in range(1, len(sys.argv)):

	n_layers = 0

	dz = 0
	dr = 0
	nz = 0
	nr = 0
	na = 0

	A_l = []
	A_z = []
	Rd_r = []
	Rd_a = []
	Tt_r = []
	Tt_a = []
	A_rz = [[]]
	Rd_ra = [[]]
	Tt_ra = [[]]

	mco = sys.argv[i]
	mco_file = open(mco)
	lines = mco_file.read().splitlines()

	stripped_lines = [l.strip() for l in lines if l.strip()!='' and l[0]!='#']

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

	if n_layers==0 or dz==0 or dr==0 or nz==0 or nr==0 or na==0:
		print('layer and detection arrays must have non-zero size ('+mco+')')
		exit()

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

	if len(A_rz) != nr:
		print('Incorrect data dimensions for A_rz ('+mco+')')
		exit()
	for l in A_rz:
		if len(l) != nz:
			print('Incorrect data dimensions for A_rz ('+mco+')')
			exit()

	if len(Rd_ra) != nr:
		print('Incorrect data dimensions for Rd_ra ('+mco+')')
		exit()
	for l in Rd_ra:
		if len(l) != na:
			print('Incorrect data dimensions for Rd_ra ('+mco+')')
			exit()

	if len(Tt_ra) != nr:
		print('Incorrect data dimensions for Tt_ra ('+mco+')')
		exit()
	for l in Tt_ra:
		if len(l) != na:
			print('Incorrect data dimensions for Tt_ra ('+mco+')')
			exit()

	folder = mco[0:-4]
	plotToPDF(folder, ['A', 'l'], A_l, [1], [''])
	plotToPDF(folder, ['A', 'z'], A_z, [dz], ['cm'])
	plotToPDF(folder, ['A', 'r', 'z'], A_rz, [dr, dz], ['cm', 'cm'])
	plotToPDF(folder, ['Rd', 'r', 'a'], Rd_ra, [dr, na/90.0], ['cm', 'deg'])
	plotToPDF(folder, ['Tt', 'r', 'a'], Tt_ra, [dr, na/90.0], ['cm', 'deg'])
	plotToPDF(folder, ['Rd', 'r'], Rd_r, [dr], ['cm'])
	plotToPDF(folder, ['Rd', 'a'], Rd_a, [na/90.0], ['deg'])
	plotToPDF(folder, ['Tt', 'r'], Tt_r, [dr], ['cm'])
	plotToPDF(folder, ['Tt', 'a'], Tt_a, [na/90.0], ['deg'])