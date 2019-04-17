#######################################################################################################
# This python script does automatic test runs of MCML.
# For specified sets of input parameters it runs all combinations.
#######################################################################################################

import sys
import os
import subprocess
import tempfile
import argparse
import datetime

from heightfieldGenerator import generateZeros, generateCurve, generateNoise



# Parameters with physical meaning from a Medical & Biological Engineering & Computing paper:
# Modelling the sampling volume for skin blood oxygenation measurements (Meglinsky, 2001)
# Table 2 Optical properties of skin model (lambda=600nm, Caucasian skin)
# k Skin layer             mus mua   g    n
# 1 stratum corneum        100 0.02  0.9  1.53
# 2 living epidermis       40  0.015 0.85 1.34
# 3 papillary dermis       30  0.07  0.8  1.4
# 4 upper blood net dermis 35  0.1   0.9  1.39
# 5 dermis                 20  0.07  0.76 1.4
# 6 deep blood net dermis  35  0.1   0.95 1.39
# 7 subcutaneous fat       15  0.03  0.8  1.44
# The refractive index of air is n=1


# Input values to test:

mua = [0.1, 1, 10] # absorption events per unit
mus = [1, 10, 100] # scattering events per unit
g = [0, 0.5, 0.9]  # scattering anisotropy https://omlc.org/classroom/ece532/class3/hg.html
n = [1.2701, 1.34, 1.44] # refractive indices
dl = [0.1, 1, 10] # layer thickness

PHOTON_MILLIONS = 100

# user can either choose one from here or use implicit flat boundaries
boundary_presets = {
	'flat50': {
		'description': 'A trivial flat boundary with 50 samples and uniform spacing 0.01.',
		'data': generateZeros(50, 0.01)
	},
	'hihill': {
		'description': 'An upward curve with high curvature',
		'data': generateCurve(0.2, 50, 0.01)
	},
	'lohill': {
		'description': 'An upward curve with low curvature',
		'data': generateCurve(0.1, 50, 0.01)
	},
	'hicrater': {
		'description': 'A downward curve with high curvature',
		'data': generateCurve(-0.2, 50, 0.01)
	},
	'locrater': {
		'description': 'A downward curve with low curvature',
		'data': generateCurve(-0.1, 50, 0.01)
	},
	'hinoise': {
		'description': 'Noise with high amplitude',
		'data': generateNoise(0.2, 50, 0.01)
	},
	'lonoise': {
		'description': 'Noise with low amplitude',
		'data': generateNoise(0.1, 50, 0.01)
	}
}


# Detection grid:

dz = 0.01 #TODO more resolution
dr = 0.01
nz = 200
nr = 500
na = 10



# Interpret options:

if len(sys.argv) < 2:
	
	print('This script tests a user-specified MCML implementation with automated parameter combinations, producing mco files and plots.')
	print('For the plots it depends on the script plotter.py, which can also be used standalone.')

	print('\nUsage: python tester.py <mcml executable> <optional executable args> -b <optional explicit boundaries(*)> -o <optional output folder> -c <optional mcml executable to compare with>')
	print('Example: python tester.py ../clustermcml-windows.exe mcmlKernel.c "-Werror" -b lohill flat50 -o mytest1')

	print('\nYou can specify an executable in any folder, which the exe will recognize as its current working directory, e.g. to find kernel files to compile.')
	
	print('\nOutput files have the following naming convention:')
	print('<number of photons>_<absorption coefficient mua>_<scattering coefficient mus>_<scattering anisotropy g>_<refractive index n>_<layer thickness dl>.mco')
	print('For each mco file, a folder with the same name is created to store the plots.')
	
	print('\n(*) Explicit boundary presets:')
	for name in boundary_presets: print(name + ' - ' + boundary_presets[name]['description'])
	
	print('\nIf no explicit boundaries are specified, implicit flat boundaries are used.')
	print('If one explicit boundary is specified, it is used above and below one layer.')
	print('If two explicit boundaries are specified, the first is used above and the second below one layer.')
	print('If more than two explicit boundaries are specified, the same layer is used multiple times so that all boundaries are used.')
	exit()



#TODO proper argument handling

# parser = argparse.ArgumentParser('Test mcml implementation with automated parameter combinations')

# parser.add_argument('mcml_exe', type=str, nargs='+', help='mcml executable')

# parser.add_argument('-b', dest='') #TODO cannot handle lists => implement new concept for providing list of explicit boundaries


mcml_exe = []
while len(sys.argv)>1 and sys.argv[1]!='-b' and sys.argv[1]!='-o' and sys.argv[1]!='-c':
	mcml_exe.append(sys.argv.pop(1))

explicit_boundaries = []
if len(sys.argv)>2 and sys.argv[1]=='-b':
	sys.argv.pop(1)
	while len(sys.argv)>1 and sys.argv[1]!='-o' and sys.argv[1]!='-c':
		b = sys.argv.pop(1)
		try:
			explicit_boundaries.append(boundary_presets[b]['data'])
		except KeyError:
			print('Preset boundary ' + b + ' does not exist')
			exit()

out_folder = os.getcwd() + '/'
if len(sys.argv)>2 and sys.argv[1]=='-o':
	sys.argv.pop(1)
	out_folder += sys.argv.pop(1) + '/'
if not os.path.exists(out_folder): os.makedirs(out_folder)

compare_mcml_exe = []
if len(sys.argv)>2 and sys.argv[1]=='-c':
	sys.argv.pop(1)
	compare_mcml_exe.append(sys.argv.pop(1))

COMPARE_MODE_ENABLED = (len(compare_mcml_exe) > 0)



# MPI:
#TODO make this optional

USE_MPI = True
mpiDebugLevel = 0
mpiMachinefile = 'mpi.txt'
mpiCmd = ['mpiexec', '-debug', str(mpiDebugLevel), '-lines', '-machinefile', mpiMachinefile, '-pwd', 'S>v2zv]08ct/lg4+9{lvi-[bcv-Qc4[|]$ne1NE~']







def makeMCIFile(mci, mco, mua_i, mus_i, g_i, n_i, dl_i, use_explicit_boundaries):
	with open(mci, 'w') as f:
		f.write('1.0\n') # file version
		f.write('1\n') # number of sims

		# sim 1
		f.write(mco+' A\n') # output filename, ASCII
		f.write(str(PHOTON_MILLIONS*1000000) + '\n') # number of photons
		f.write(str(dz)+' '+str(dr)+'\n') # detection deltas in cm (dz, dr)
		f.write(str(nz)+' '+str(nr)+' '+str(na)+'\n') # detection array dims (z, r, a)

		# number of layers
		if use_explicit_boundaries and len(explicit_boundaries) > 1:
			f.write(str(len(explicit_boundaries)-1) + '\n')
		else:
			f.write('1\n')

		f.write('1.0\n') # n above

		if use_explicit_boundaries and len(explicit_boundaries) > 1:
			# In the current "boundary-centric" test, one layer is repeated
			# with the same parameters to use all specified boundaries
			for i in range(0, len(explicit_boundaries)-1):
				boundary_indicator = 'b ' + str(int(len(explicit_boundaries[i]) / 2)) + ' '
				f.write(boundary_indicator + ' '.join(str(f) for f in explicit_boundaries[i]) + '\n') # boundary
				f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl_i)+'\n') # layer
			boundary_indicator = 'b ' + str(int(len(explicit_boundaries[len(explicit_boundaries)-1]) / 2)) + ' '
			f.write(boundary_indicator + ' '.join(str(f) for f in explicit_boundaries[len(explicit_boundaries)-1]) + '\n') # boundary
		elif use_explicit_boundaries and len(explicit_boundaries) > 0:
			boundary_indicator = 'b ' + str(int(len(explicit_boundaries[0]) / 2)) + ' '
			f.write(boundary_indicator + ' '.join(str(f) for f in explicit_boundaries[0]) + '\n') # boundary
			f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl_i)+'\n') # layer
			f.write(boundary_indicator + ' '.join(str(f) for f in explicit_boundaries[0]) + '\n') # boundary
		else:
			f.write(str(n_i)+' '+str(mua_i)+' '+str(mus_i)+' '+str(g_i)+' '+str(dl_i)+'\n') # layer 1

		f.write('1.0') # n below



mcos = [] # gets filled with names of produced output files
failed_mcos = []

simcount = 0

def runOneTest(mua_i, mus_i, g_i, n_i, dl_i, maxsims):

	global simcount
	simcount += 1

	# Make a compact string containing important mcml config paramenters of current run
	config_str = out_folder + str(PHOTON_MILLIONS)+'M_'+str(mua_i)+'mua_'+str(mus_i)+'mus_'+str(g_i)+'g_'+str(n_i)+'n_'+str(dl_i)+'dl'

	# Compare mode requires 2 runs, for both executables
	C = (2 if COMPARE_MODE_ENABLED else 1)
	is_compare = 1 # index of the comparison run

	for c_i in range(0, C):

		mco = config_str

		# In compare mode append executable basenames to mco for clear differentiation
		if c_i==is_compare:
			mco += '_' + os.path.splitext(os.path.basename(compare_mcml_exe[0]))[0]
		elif COMPARE_MODE_ENABLED:
			mco += '_' + os.path.splitext(os.path.basename(mcml_exe[0]))[0]

		mco = os.path.abspath(mco + '.mco')

		if (os.path.isfile(mco)):
			print('[SKIP] existing ' + str(mco))
			mcos.append(mco) # need this later for plotting
			continue

		mci = os.path.join(tempfile.mkdtemp(), 'test.mci')

		use_explicit_boundaries = (len(explicit_boundaries) > 0 and c_i!=is_compare) # never want explicit boundaries for compare exe
		makeMCIFile(mci, mco, mua_i, mus_i, g_i, n_i, dl_i, use_explicit_boundaries)

		print('===========================================================')
		print('[STARTING] Simulation ' + str(simcount) + ' of ' + str(maxsims) +
			', compare = ' + str(c_i==is_compare) +
			', t = ' + str(datetime.datetime.now()) )
		print('===========================================================')
		print('[INPUT] ' + mci + ':')
		with open(mci, 'r') as f:
			print(f.read())
		print('===========================================================')

		success = False
		retry_count = 0
		max_retries = 3

		while not success and not retry_count > max_retries:

			# If using MPI, we actually execute mpiexec and pass mcml executable as argument

			# First part of executable path is used as cwd, second as mpiexec argument
			pathPair = os.path.split(mcml_exe[0])

			cwd = pathPair[0]

			cmd = [*mcml_exe, mci]

			if USE_MPI:
				cmd[0] = pathPair[1]
				cmd = mpiCmd + cmd

			# Compare exe currently never uses MPI or another cwd
			if c_i==is_compare:
				cmd = [*compare_mcml_exe, mci]
				cwd = '.'

			print('[EXEC] ' + str(cmd) + ':')
			result = subprocess.run(cmd, cwd=cwd)

			print('===========================================================')
			try:
				result.check_returncode()

				print('[OUTPUT] ' + mco)


				if COMPARE_MODE_ENABLED and c_i==is_compare:
					print('[PLOTTING]', mcos[len(mcos)-1], mco)
					subprocess.run(['python', 'plotter.py', mcos[len(mcos)-1], mco, '--compare', '-o', config_str])
				elif not COMPARE_MODE_ENABLED:
					print('[PLOTTING]')
					subprocess.run(['python', 'plotter.py', mco])

				mcos.append(mco)

				success = True

			except subprocess.CalledProcessError as err:
				print('[FAILED] ' + str(err))
				print('Retrying...')
				retry_count += 1
				if retry_count == max_retries:
					failed_mcos.append(mco)
			print('\n\n')



# Uncomment to test all combinations of above parameters
# for mua_i in mua:
# 	for mus_i in mus:
# 		for g_i in g:
# 			for n_i in n:
# 				for dl_i in dl:

# 					runOneTest(mua_i, mus_i, g_i, n_i, dl_i, len(mua)*len(mus)*len(g)*len(n)*len(dl))


maxsims = 13
# Varying mua, mus
runOneTest(1, 1, 0, 1.3, 0.1, maxsims)
runOneTest(1, 10, 0, 1.3, 0.1, maxsims)
runOneTest(10, 100, 0, 1.3, 0.1, maxsims)
# Varying g
runOneTest(1, 10, -0.5, 1.3, 0.1, maxsims)
runOneTest(1, 10, 0, 1.3, 0.1, maxsims)
runOneTest(1, 10, 0.5, 1.3, 0.1, maxsims)
# Varying n
runOneTest(1, 10, 0, 1.1, 0.1, maxsims)
runOneTest(1, 10, 0, 1.3, 0.1, maxsims)
runOneTest(1, 10, 0, 1.5, 0.1, maxsims)
runOneTest(1, 10, 0, 2.0, 0.1, maxsims)
# Varying d
runOneTest(1, 10, 0, 1.3, 0.01, maxsims)
runOneTest(1, 10, 0, 1.3, 0.1, maxsims)
runOneTest(1, 10, 0, 1.3, 1, maxsims)




print('[FINISHED] ' + str(len(mcos)) + ' mco files successfully produced.')
print(str(len(failed_mcos)) + ' simulations failed' + (':' if len(failed_mcos) > 0 else '.'))
for mco in failed_mcos: print(mco)